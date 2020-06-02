package upv.es.bidding;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import genius.core.Bid;
import genius.core.BidHistory;
import genius.core.Domain;
import genius.core.bidding.BidDetails;
import genius.core.boaframework.BOAparameter;
import genius.core.boaframework.NegotiationSession;
import genius.core.boaframework.OpponentModel;
import genius.core.boaframework.SortedOutcomeSpace;
import genius.core.issue.ISSUETYPE;
import genius.core.issue.Issue;
import genius.core.issue.IssueDiscrete;
import genius.core.issue.Value;
import genius.core.issue.ValueDiscrete;
import genius.core.issue.ValueInteger;
import genius.core.issue.ValueReal;
import genius.core.utility.AbstractUtilitySpace;
import genius.core.utility.AdditiveUtilitySpace;
import negotiator.boaframework.opponentmodel.tools.UtilitySpaceAdapter;


public class PasqualOpponentModel extends OpponentModel {
	
	protected int verbose = 1;
    
    // Global frequency model
	protected Hashtable<String, Object> freqModel;
	protected Hashtable<Number, String> issueIdToName;
    protected Hashtable<String, Double> issueWeight;
    // Bid windows
	protected List<Bid> win1;
	protected List<Bid> win2;

	protected SortedOutcomeSpace sortedOutcomeSpace;

	protected int WINDOW_SIZE = 5;
	protected double THRESHOLD = 0.1;
	protected double ALPHA = 10.0;
	protected double BETA = 5.0;
	
	
	@Override
	public String getName() {
		return "Pasqual Frequency Model";
	}
	
	/** ====================================================================
	 * 	========================== INITIALIZATION ==========================
	 *  ====================================================================
	 */
	
	@Override
	// Opponent Modeling initialization
	public void init(NegotiationSession negotiationSession, Map<String, Double> parameters) {
		super.init(negotiationSession, parameters);
		
		// Initialize window
		win1 = new ArrayList<Bid>();
		win2 = new ArrayList<Bid>();
		if(parameters.containsKey("window_size")) {
			WINDOW_SIZE = (int) parameters.get("window_size").intValue();
		}
		if(parameters.containsKey("threshold")) {
			THRESHOLD = (double) parameters.get("threshold");
		}
		if(parameters.containsKey("alpha")) {
			ALPHA = (double) parameters.get("alpha");
		}
		if(parameters.containsKey("beta")) {
			BETA = (double) parameters.get("beta");
		}
		
		/** 
		 * Initialize Frequency Model 
		 **/
		if (verbose > 0) {
			System.out.println("====================================================================");
			System.out.println("========================== INITIALIZATION ==========================");
			System.out.println("====================================================================");
			System.out.println();
			System.out.printf("Verbosity level %d\n\n", verbose);
			System.out.println("-------------------------------\nInitializing frequency model");
		}
		
		// Initialize tables
		freqModel = new Hashtable<String, Object>();
		issueIdToName = new Hashtable<Number, String>();
		issueWeight = new Hashtable<String, Double>();
		// Get list of domain issues
		List<Issue> issues = negotiationSession.getIssues();
		double initialWeight = 1.0/issues.size();
		
		for (Issue i : issues) {
			
			// Map issue number to its name
			issueIdToName.put(i.getNumber(), i.getName());
			
			// Map issue name to its frequency model weight
			issueWeight.put(i.getName(), initialWeight);
			
			if (verbose > 1) System.out.println("Issue: "+i.getName()+" is of type "+i.getType());
			
			if (i.getType() == ISSUETYPE.DISCRETE) {
				
				IssueDiscrete d = (IssueDiscrete) i;
				
				// Initialize appearances of every possible value of the issue to 1
				Hashtable<String,Integer> countValues = getAllValues(d, 1);
				freqModel.put(i.getName(), countValues);
			}
			else if (i.getType() == ISSUETYPE.REAL) {
				freqModel.put(i.getName(), new ArrayList<Double>());
			} 
			else if (i.getType() == ISSUETYPE.INTEGER) {
				freqModel.put(i.getName(), new ArrayList<Integer>());
			}			
			
		}
		if (verbose > 0) {
			System.out.println("Frequency model initialized: "+freqModel.toString());
			System.out.println("Issues id to name mapping: "+issueIdToName.toString());
			System.out.println("Initial weights by issue: "+printWeights()+"\n-------------------------------");
			System.out.println();
		}
		
	}
	
	/**
	 * Initializes the frequency model for a given issue using the Laplace Smoothing.
	 * Every possible value for that issue will be assigned a 1, as if it had already
	 * been seen in an offer once.
	 */
	public Hashtable<String,Integer> getAllValues(IssueDiscrete i, int alpha) {
		// Dictionary to map every value with the amount of times it has appeared
		// in an opponent's bid
		Hashtable<String,Integer> countValues = new Hashtable<String,Integer> ();
		
		// List of possible values for the issue
		List<ValueDiscrete> possibleValues = i.getValues();
		
		for (ValueDiscrete v : possibleValues) {
			countValues.put(v.getValue(), alpha);
		}
		
		return countValues;
	}
	
	/** ====================================================================
	 * 	========================= FREQUENCY MODEL ==========================
	 *  ====================================================================
	 */
	
	
	/**
	 * Given a bid, increases the counters of frequency model
	 * for the issue values contained in such bid
	 * @param bid
	 */
	@SuppressWarnings("unchecked")
	protected void updateFrequencyModel(Bid bid) {
		HashMap<Integer, Value> values = bid.getValues();
		if (verbose > 0) {
			System.out.println("-------------------------------\nUpdating frequency model...");
			System.out.println("Bid values: "+values.toString());
		}
		Set<Number> keys = issueIdToName.keySet();
		// for every issue
		for (Number key : keys) {
			// get issue name and its value
			String name = issueIdToName.get(key);
			Value value = values.get(key);
			
			if (verbose > 1) System.out.println("Issue "+name+" has value "+value);
			
			// Get object that the frequency model has stored for the corresponding issue
			// Depending on issue type it can be a Hashtable or an ArrayList of numbers
			Object o = freqModel.get(name);
			
			// If the issue is a Discrete issue
			if (o instanceof Hashtable) {
				// Update number of appearances of its value
				Hashtable<String, Integer> countValues = (Hashtable<String, Integer>) o ; 
				if (countValues.containsKey(value.toString())) {
					int count = countValues.get(value.toString());
					count++;
					// update hashtable
					countValues.replace(value.toString(), count);
				}
				else {
					System.out.printf("ERROR :: Frequency model found value {%s} for issue {%s} that wasn't initialized\n", name, value.toString());
				}
				// update hashtable
				freqModel.replace(name, countValues);
			}
			// If the issue is an Integer issue
			else if (value.getType() == ISSUETYPE.INTEGER) {
				System.out.println("WARNING :: Frequency model not supposed for non-discrete issues -------------------");
				// Add value to list of stored values
				ArrayList<Integer> valueList = (ArrayList<Integer>) o ; 
				ValueInteger iv = (ValueInteger) value;
				// valueList.add(iv.getValue());
			}
			// If the issue is a Real issue
			else if (value.getType() == ISSUETYPE.REAL) {
				System.out.println("WARNING :: Frequency model not supposed for non-discrete issues -------------------");
				// Add value to list of stored values
				ArrayList<Double> valueList = (ArrayList<Double>) o ; 
				ValueReal iv = (ValueReal) value;
				//valueList.add(iv.getValue());
			}
			
		} // end of keys For
		if (verbose > 0) System.out.println("Updated freq. model: "+freqModel.toString()+"\n-------------------------------\n");
		System.out.println("++++++++++++The opponent has sent "+negotiationSession.getOpponentBidHistory().getHistory().size()+ " offers so far++++++++++");
		
	}
	
	/**
	 * ****************************** CURRENTLY NOT IN USE *******************************
	 * Compares the global frequency of every issue value against the frequency in the window.
	 * If the frequencies differ, we will assume that the opponent is conceding on that issue,
	 * and update the weights of the issues in which it does not concede
	 * ****************************** CURRENTLY NOT IN USE *******************************
	 */
	protected void updateWeightsGlobal() {
		if (verbose > 0) {
			System.out.println("-------------------------------\nUpdating model weights");
			System.out.println("Current weights: "+issueWeight.toString());
		}
		int totalOffers = negotiationSession.getOpponentBidHistory().getHistory().size();
		Hashtable<String, Hashtable<String, Double>> globalDic = getFrequency(freqModel, totalOffers);
		Hashtable<String, Hashtable<String, Double>> windowDic = getWindowFreq(win1);
		Hashtable<String, Boolean> concealed = new Hashtable<String, Boolean>();
		
		if (verbose > 0) {
			System.out.println("Global dic "+globalDic.toString());
			System.out.println("Window dic "+windowDic.toString());
		}
		
		Set<String> issues = globalDic.keySet();
		for (String issue : issues) {
			
			if (verbose > 0) {
				System.out.println("Checking frequencies for issue "+issue);
			}
			
			// Get global and window frequencies for every value
			
			Hashtable<String, Double> globalValues = globalDic.get(issue);
			Hashtable<String, Double> windowValues = windowDic.get(issue);
			
			if (verbose > 1) {
				System.out.println("{\n\tGlobal frequencies "+globalValues.toString());
				System.out.println("\tWindow frequencies "+windowValues.toString()+"\n}");
			}
			
			Set<String> values = globalValues.keySet();
			
			double maxGlobalFreq = -1;
			String maxGlobalValue = "x";
			
			double maxWindowFreq = -1;
			String maxWindowValue = "y";
			
			// Get the issue value with maximum frequency globaly and inside the window
			for (String value : values) {
				
				double globalFreq = globalValues.get(value);
				if (globalFreq > maxGlobalFreq) {
					maxGlobalFreq = globalFreq;
					maxGlobalValue = value;
				}
				
				double windowFreq = windowValues.get(value);
				if (windowFreq > maxWindowFreq) {
					maxWindowFreq = windowFreq;
					maxWindowValue = value;
				}
				
			}
			
			// If the value with maximum frequency changed, 
			// the opponent has conceded on the issue
			if (!maxGlobalValue.equals(maxWindowValue)) {
				if (verbose > 1) {
					System.out.printf("The maximum values do not match : %s != %s\n", maxGlobalValue, maxWindowValue);
				}
				concealed.put(issue, true);
			} 
			// If it did not change but its difference is grater than a given threshold,
			// we assume that the opponent has conceded on the issue
			else if (Math.abs(maxGlobalFreq - maxWindowFreq) > THRESHOLD) {
				if (verbose > 1) {
					System.out.printf("The difference between maximum values overcomes threshold : | %f - %f | > %f\n", 
							maxGlobalFreq, maxWindowFreq, THRESHOLD);
				}
				concealed.put(issue, true);
			}
			else {
				if (verbose > 1) {
					System.out.printf("The maximum frequencies match : %s == %s, | %f - %f | < %f\n", 
							maxGlobalValue, maxWindowValue, maxGlobalFreq, maxWindowFreq, THRESHOLD);
				}
				concealed.put(issue, false);
			}
			
		} // end of issues For
		
		// For every issue in which opponent has not conceded, increase its weight
		// according to current negotiation time
		if (verbose > 1) {
			System.out.println("Concealed dic: "+concealed.toString());
		}
		for (String issue : issues) {
			if (!concealed.get(issue)) {
				double oldW = issueWeight.get(issue);
				double newW = oldW + getDelta(ALPHA, BETA);
				issueWeight.replace(issue, newW);
				if (verbose > 1) System.out.printf("Increasing weight of issue %s from %f to %f\n", issue, oldW, newW);
			}
		}
		// Rescale issue weights so that they add 1
		scaleWeights();
		
		if (verbose > 0) {
			System.out.println("Updated weights: "+issueWeight.toString()+"\n-------------------------------\n");
		}
	}
	
	/**
	 * Compares old and new windows and checks, for every issue, if the opponent has concealed.
	 * This is done by comparing the change in the frequencies for every possible value of an issue
	 * between both windows.
	 * If the value with maximum frequency for one issue does not match for both windows, we 
	 * consider that the opponent concealed.
	 * The weigths for the non-concealed issues are increased at the end of the comparison
	 */
	protected void updateWeights() {
		if (verbose > 0) {
			System.out.println("-------------------------------\nUpdating model weights");
			System.out.println("Current weights: "+issueWeight.toString());
		}
		int totalOffers = negotiationSession.getOpponentBidHistory().getHistory().size();
		Hashtable<String, Hashtable<String, Double>> win1Dic = getWindowFreq(win1);
		Hashtable<String, Hashtable<String, Double>> win2Dic = getWindowFreq(win2);
		Hashtable<String, Boolean> concealed = new Hashtable<String, Boolean>();
		
		if (verbose > 0) {
			System.out.println("Window 1 dic "+win1Dic.toString());
			System.out.println("Window 2 dic "+win2Dic.toString());
		}
		
		Set<String> issues = win1Dic.keySet();
		for (String issue : issues) {
			
			if (verbose > 0) {
				System.out.println("Checking frequencies for issue "+issue);
			}
			
			// Get global and window frequencies for every value
			
			Hashtable<String, Double> win1Values = win1Dic.get(issue);
			Hashtable<String, Double> win2Values = win2Dic.get(issue);
			
			if (verbose > 1) {
				System.out.println("{\n\tWindow 1 frequencies "+win1Values.toString());
				System.out.println("\tWindow 2 frequencies "+win2Values.toString()+"\n}");
			}
			
			Set<String> values = win1Values.keySet();
			
			double maxWin1Freq = -1;
			String maxWin1Value = "x";
			
			double maxWin2Freq = -1;
			String maxWin2Value = "y";
			
			// Get the issue value with maximum frequency globaly and inside the window
			for (String value : values) {
				
				double win1Freq = win1Values.get(value);
				if (win1Freq > maxWin1Freq) {
					maxWin1Freq = win1Freq;
					maxWin1Value = value;
				}
				
				double win2Freq =win2Values.get(value);
				if(win2Freq > maxWin2Freq) {
					maxWin2Freq = win2Freq;
					maxWin2Value = value;
				}
				
			}
			
			// If the value with maximum frequency changed, 
			// the opponent has conceded on the issue
			if (!maxWin2Value.equals(maxWin1Value)) {
				if (verbose > 1) {
					System.out.printf("The maximum values do not match : %s != %s\n", maxWin1Value, maxWin2Value);
				}
				concealed.put(issue, true);
			} 
			// If it did not change but its difference is grater than a given threshold,
			// we assume that the opponent has conceded on the issue
			else if (Math.abs(maxWin1Freq - maxWin2Freq) > THRESHOLD) {
				if (verbose > 1) {
					System.out.printf("The difference between maximum values overcomes threshold : | %f - %f | > %f\n", 
							maxWin1Freq, maxWin2Freq, THRESHOLD);
				}
				concealed.put(issue, true);
			}
			else {
				if (verbose > 1) {
					System.out.printf("The maximum frequencies match : %s == %s, | %f - %f | < %f\n", 
							maxWin1Value, maxWin2Value, maxWin1Freq, maxWin2Freq, THRESHOLD);
				}
				concealed.put(issue, false);
			}
			
		} // end of issues For
		
		// For every issue in which opponent has not conceded, increase its weight
		// according to current negotiation time
		if (verbose > 1) {
			System.out.println("\nConcealed dic: "+concealed.toString()+"\n");
		}
		for (String issue : issues) {
			if (!concealed.get(issue)) {
				double oldW = issueWeight.get(issue);
				double newW = oldW + getDelta(10, 5);
				issueWeight.replace(issue, newW);
				if (verbose > 1) System.out.printf("Increasing weight of issue %s from %f to %f\n", issue, oldW, newW);
			}
		}
		// Rescale issue weights so that they add 1
		scaleWeights();
		
		if (verbose > 0) {
			System.out.println("Updated weights: "+issueWeight.toString()+"\n-------------------------------\n");
		}
		
	}
	
	
	/**
	 * Updates frequency model and adds the received opponent bid
	 * to the corresponding window
	 */
	@Override
	protected void updateModel(Bid bid, double time) {
		
		updateFrequencyModel(bid);		
		
		
		if(win2.size() == WINDOW_SIZE) {
			if (verbose > 0) System.out.println("Window reached size "+WINDOW_SIZE);
			
			// Update weights
			updateWeights();
			
			// Copy win2 into win1
			win1.clear();
			for (Bid b : win2) {
				Bid copy = new Bid(b);
				win1.add(copy);
			}
			
			// Clear win2
			win2.clear();		
		}
		
		// Add bids to window 1 (only at the beginning of the negotiation
		if(win1.size() < WINDOW_SIZE) {
			win1.add(bid);
		}
		// If window 1 already full, add bids to window 2
		else if(win2.size() < WINDOW_SIZE) {
			win2.add(bid);
		}
		
		System.out.printf("Round %.0f weights: %s", negotiationSession.getTimeline().getCurrentTime(), printWeights());
	}
	
	/**
	 * Calculates the opponent's utility of a given bid
	 * Taking into account the global frequency model
	 * @param bid
	 * @return
	 */
	@SuppressWarnings("unchecked")
	protected double getOpponentUtility(Bid bid) {
		HashMap<Integer, Value> values = bid.getValues();
		
		//double [] estimation = new double[values.size()];
		double utility = 0;
		
		
		Set<Number> keys = issueIdToName.keySet();
		// for every issue
		for (Number key : keys) {
			
			// get issue name and its value
			String name = issueIdToName.get(key);
			Value value = values.get(key);
			double count = -1;
			double max = -1;
			double temp;
			
			Hashtable<String,Integer> countValues = (Hashtable<String,Integer>) freqModel.get(name);
			
			// Variable "count" stores the number of times value 
			// "value" has appeared in a bid since the beginning
			count = (double) countValues.get(value.toString());
			
			// Variable "max" stores the number of times that the maximum value 
			// for that issue appeared in a bid since the beginning
			Set<String> seenValues = countValues.keySet();
			for (String key1 : seenValues) {
				temp = countValues.get(key1);
				if (temp > max)	max = temp;
			}
			
			double valueFreq = count * 1.0 / max;
			double issueScore = issueWeight.get(name) * valueFreq;
			utility += issueScore;

			if (verbose > 1) {
				System.out.println("{");
				System.out.println("\tValue "+value.toString()+" appears "+count+" times in issue "+key);
				System.out.println("\tMax times a value for issue "+key+" appeared is "+max);
				System.out.println("\tThe weight for issue "+key+" is "+issueWeight.get(key));
				System.out.println("\tThe estimate for "+value.toString()+" is "+issueScore);
				System.out.println("}");
			}
			
		} // end of issues For
		
		return utility;
	}
	
	/**
	 * Returns an estimation of the opponent's utility for a given bid.
	 * Initially returns the agent's own utility.
	 * Once there are enough bids in the frequency model, the estimation 
	 * will be based on it. 
	 */
	@Override
	// Estimation of the opponent's utility of a given bid
	public double getBidEvaluation(Bid bid) {
		if (verbose > 1) System.out.println("-------------------------------\nEvaluating bid "+bid.getValues().toString());
		// CHANGE
		
		double utility = 0;
		
		if (negotiationSession.getTimeline().getCurrentTime() > 5) {
		
			utility = getOpponentUtility(bid);
			
			if (verbose > 1) {
				System.out.printf("Utility: %3.5f\n-------------------------------\n\n", utility);
			}
			
			return utility;
		} 
		else {
			if (verbose > 1) {
				System.out.println("Not enough opponent bids, returning own utility\n-------------------------------");
			}
			return negotiationSession.getUtilitySpace().getUtility(bid);
		}
		
	}
	
	@Override
	public AdditiveUtilitySpace getOpponentUtilitySpace() {	
		AdditiveUtilitySpace utilitySpace = new UtilitySpaceAdapter(this, this.negotiationSession.getDomain());
		return utilitySpace;
	}
	
	@Override
	// Parameters for the Opponent Model component in Genius
	public Set<BOAparameter> getParameterSpec(){
		Set<BOAparameter> set = new HashSet<BOAparameter>();
		set.add(new BOAparameter("window_size", 5.0, "Window capacity in number of opponent's bid"));
		set.add(new BOAparameter("threshold", 0.1, "Min. difference in value frequency to consider the opponent concealed in an issue"));
		set.add(new BOAparameter("alpha", 10.0, "Base weight increase"));
		set.add(new BOAparameter("beta", 5.0, "Controls influence of time in weight update. Higher values = more increment"));
		return set;
	}


	/** ====================================================================
	 * 	======================= AUXILIARY FUNCTIONS ========================
	 *  ====================================================================
	 */
	
	/**
	 * Given a hash table with value appearances, returns a similar structure 
	 * with value frequencies with respect to a specified total
	 */
	public Hashtable<String,Hashtable<String,Double>> getFrequency(Hashtable<String, Object> model, int total) {
		
		// Hashtable to store the frequency of every value of every issue
		Hashtable<String, Hashtable<String,Double>> frequencies = new Hashtable<String, Hashtable<String,Double>>();
		
		Set<String> keys = model.keySet();
		for (String issue : keys) {
			
			// Hashtable to store the frequency of every value of a specific issue
			Hashtable<String,Double> valueFreq = new Hashtable<String,Double>();
			
			Object o = model.get(issue);
			// If the issue is a Discrete issue
			if (o instanceof Hashtable) {
				@SuppressWarnings("unchecked")
				Hashtable<String,Integer> countValues = (Hashtable<String,Integer>) o;		
			
				Set<String> possibleValues = countValues.keySet();
				for (String value : possibleValues) {
					// Get number of times the value appeared in a bid
					double count = countValues.get(value);
					// Frequency is calculated as "count" divided by a given total (number of rounds, elements of a window)
					double freq = count * 1.0 / (total+possibleValues.size());
					// Store freq in hashtable
					valueFreq.put(value, freq);
				}
				frequencies.put(issue, valueFreq);
			}
		}
		
		return frequencies;
	}
	
	/**
	 * Given a window of bids, returns a Hash table with the frequency
	 * of every issue value inside that window
	 */
	public Hashtable<String, Hashtable<String,Double>> getWindowFreq(List<Bid> window) {
		
		// Get base empty model with Laplace smoothing
		Hashtable<String, Object> model = getEmptyModel();
		
		// For every bid in the window
		for (Bid bid : window) {
			HashMap<Integer, Value> bidValues = bid.getValues();
			
			Set<Number> keys = issueIdToName.keySet();
			// For every issue
			for (Number key : keys) {
				
				// Get issue name and its value
				String name = issueIdToName.get(key);
				Value value = bidValues.get(key);

				
				// Add one to its appearances
				Hashtable<String,Integer> countValues = (Hashtable<String,Integer>) model.get(name);
				int count = countValues.get(value.toString());
				count += 1;
				countValues.replace(value.toString(), count);
			}
		}
		
		// Once the model is complete with the window data, calculate frequencies and return
		return getFrequency(model, window.size());		
	}
	
	/**
	 * Returns an empty dictionary containing, for every issue, a dictionary
	 * with all its possible values and their appearances, which are set 
	 * to 1 so as to apply Laplace Smoothing
	 * @return
	 */
	public Hashtable<String, Object> getEmptyModel() {
		// Initialize tables
		Hashtable<String, Object> emptyModel = new Hashtable<String, Object>();

		// Get list of domain issues
		List<Issue> issues = negotiationSession.getIssues();
		double initialWeight = 1.0/issues.size();
		
		for (Issue i : issues) {
			String name = i.getName();				
			IssueDiscrete d = (IssueDiscrete) i;				
			// Initialize appearances of every possible value of the issue to 1
			Hashtable<String,Integer> countValues = getAllValues(d, 1);
			emptyModel.put(i.getName(), countValues);	
			
		}
		
		return emptyModel;
	}
	
	/**
	 * Returns the amount in which to increase an issue weight when updated.
	 * It takes into account the negotiation time.
	 * @param alpha
	 * @param beta
	 * @return
	 */
	protected double getDelta(double alpha, double beta) {
		return alpha * (1 - Math.pow(negotiationSession.getTime(), beta));
	}
	
	/**
	 * Scales the weights of the issues so that they add up to 1
	 */
	protected void scaleWeights() {
		Set<String> issues = issueWeight.keySet();
		double sum = 0;
		for (String issue : issues) {
			sum += issueWeight.get(issue);
		}
		for (String issue : issues) {
			double w = issueWeight.get(issue);
			w /= sum;
			issueWeight.replace(issue, w);
		}
	}
	
	/**
	 * Helper method to add all elements in an array
	 */
	 public double sum(double...values) {
	   double result = 0;
	   for (double value:values)
	     result += value;
	   return result;
	}
	 
	/**
	 * Helper method to visualize an array with its values
	 * in the standard output
	 */
	public String printArray(double [] a) {
		String s = "[ ";
		for (int i = 0; i < a.length; i++) {
			if (i == a.length-1) {
				s += a[i]+" ]";
			} else {
				s += a[i]+", ";
			}
		}
		return s;
	}
	
	/**
	 * Helper method to visualize the issue weight dictionary
	 * in the standard output
	 * @return
	 */
	protected String printWeights() {
		String s = "{\n";
		Set<String> issues = issueWeight.keySet();
		for (String issue : issues) {
			String issueName = String.format("%30s", issue);
			String rounded = String.format("%.8f", issueWeight.get(issue));
			
			s += "\t"+issueName+" : "+rounded+"\n";
		}
		s += "}\n";		
		return s;
	}
	
}
