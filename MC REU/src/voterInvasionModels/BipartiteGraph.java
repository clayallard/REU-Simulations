package voterInvasionModels;

import static org.junit.jupiter.api.DynamicTest.stream;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import java.util.Random;

/**
 * The bipartite graph is a graph that consists of two groups. One group of size
 * n and the other of size m. Every member of the n sized group has a connection
 * to every member of the m sized group and vice versa.
 * 
 * In this, we simulate the voter and invasion models on the complete bipartite
 * graph. Each person in a population is given opinion 0 or 1. We randomly
 * select a person in the population (person A), then randomly select one of
 * their neighbors (person B). In the voter model, person A adopts the opinion
 * of person B. In the invasion model, person A forces their opinion onto person
 * B. At each time step, we observe the number of 1 opinions of the large group
 * N and the number of 1 opinions of the small group M. This goes on until there
 * is a consensus in the population. The purpose of this project is to make
 * observations such as the average time to consensus, or the quasi-stationary
 * distribution of the process. That is, conditioned on not being at a
 * consensus, the distribution of the amount of 1 opinions in each group.
 * 
 * @author Clayton Allard
 *
 */
public class BipartiteGraph {

	// list of people.
	private ArrayList<Person> people = new ArrayList<>();
	// determines if voter or invasion model.
	private boolean voterModel;
	// current amount of time steps.
	private int time;
	// amount of people in the small group M.
	private int smallGroupAmount;
	// amount of people in the large group N.
	private int largeGroupAmount;
	// current amount of people with opinion 1 in small group and large group (i.e.
	// (4, 3)).
	private int[] currState = new int[2];
	// list of states from each time step (Ex:[(2,1), (2,2), (2,2), (1,2)]).
	private ArrayList<String> stateList = new ArrayList<>();

	/**
	 * Constructs graph consisting of two groups. A large group and a small group.
	 * Everyone in the large group is connected to everyone in the small group and
	 * vice versa. But no one in the same group is connected to eachother.
	 * 
	 * @param largeGroup - number of people in the large group.
	 * @param smallGroup - number of people in the large group.
	 * @param voter      - voter model if true. Invasion model if false.
	 */
	public BipartiteGraph(int largeGroup, int smallGroup, boolean voter) {
		voterModel = voter;
		smallGroupAmount = Math.min(smallGroup, largeGroup);
		largeGroupAmount = Math.max(smallGroup, largeGroup);
		Random rng = new Random();

		// add people to the graph
		for (int i = 1; i <= largeGroup + smallGroup; i++) {
			people.add(new Person(i, rng.nextInt(2)));
		}

		// keep track of amount of 1 in small group.
		for (int i = 0; i < smallGroupAmount; i++) {
			if (people.get(i).getOpinion() == 1) {
				currState[1]++;
			}
		}

		// keep track of amount of 1 in large group.
		for (int j = smallGroupAmount; j < largeGroup + smallGroup; j++) {
			if (people.get(j).getOpinion() == 1) {
				currState[0]++;
			}
		}

		// set appropriate relations
		for (int i = 0; i < smallGroupAmount; i++) {
			for (int j = smallGroupAmount; j < largeGroup + smallGroup; j++) {
				addRelation(people.get(i), people.get(j));
			}
		}
		people.sort((x, y) -> x.getID() - y.getID());
	}

	/**
	 * Adds connection between two people.
	 * 
	 * @param p1 - person 1.
	 * @param p2 - person 2.
	 */
	private void addRelation(Person p1, Person p2) {
		if (p1.equals(p2)) {
			return;
		}
		p1.addNeighbor(p2);
	}

	/**
	 * Resets the graph to its default state of starting at time 0.
	 */
	private void reset() {
		time = 0;
		stateList = new ArrayList<>();
	}

	/**
	 * Manually set the current state of the graph. That is, set the amount of
	 * people with opinion 1 in the large group and small group.
	 * 
	 * @param k
	 * @param h
	 */
	public void set(int k, int h) {
		reset();
		currState[0] = k;
		currState[1] = h;
		// set opinions for small group.
		for (int i = 0; i < smallGroupAmount; i++) {
			if (i < h) {
				people.get(i).setOpinion(1);
			} else {
				people.get(i).setOpinion(0);
			}
		}
		// set opinions for small group.
		for (int i = smallGroupAmount; i < people.size(); i++) {
			if (i < k + smallGroupAmount) {
				people.get(i).setOpinion(1);
			} else {
				people.get(i).setOpinion(0);
			}
		}
//		stateList.add((double) currState[0] / largeGroupAmount + "\t" + (double) currState[1] / smallGroupAmount);
	}

	/**
	 * Takes step in the process.
	 * 
	 * @return string representing who took over who's opinion.
	 */
	public String step() {
		Random rng = new Random();
		return step(rng);
	}

	/**
	 * Takes step in the process.
	 * 
	 * @param rng - to set seed.
	 * @return string representing who took over who's opinion.
	 */
	public String step(Random rng) {
		time++;
		if (voterModel) {
			return voterStep(rng);
		} else {
			return invasionStep(rng);
		}
	}

	/**
	 * Run the process until there is a consensus.
	 */
	public void simulation() {
		simulation(new Random());
	}

	/**
	 * Run the process until there is a consensus.
	 * 
	 * @param rng - to set seed.
	 */
	public void simulation(Random rng) {
		while (!consensus())
			step(rng);
	}

	/**
	 * Amount of steps the process has been running.
	 */
	public int time() {
		return time;
	}

	/**
	 * Taking a step with the voter model. We randomly select a person who takes on
	 * the opinion of one of their randomly selected neighbors.
	 * 
	 * @param rng - to set seed
	 * @return string representing who took over who's opinion.
	 */
	private String voterStep(Random rng) {
		// get the correct people.
		Person changeOpinion = people.get(rng.nextInt(people.size()));
		ArrayList<Person> neighbors = changeOpinion.neighbors();
		Person forceOpinion = neighbors.get(rng.nextInt(neighbors.size()));
		// see if the forced opinion comes from small group.
		if (forceOpinion.getID() <= smallGroupAmount) {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[0]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[0]--;
			}
			// see if forced opinion comes from large group.
		} else {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[1]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[1]--;
			}
		}

		changeOpinion.setOpinion(forceOpinion.getOpinion());
//		stateList.add((double) currState[0] / largeGroupAmount + "\t" + (double) currState[1] / smallGroupAmount);
		return changeOpinion.getID() + " <- " + forceOpinion.getID();
	}

	/**
	 * Taking a step with the voter model. We randomly select a person who takes on
	 * the opinion of one of their randomly selected neighbors.
	 * 
	 * @param rng - to set seed
	 * @return string representing who took over who's opinion.
	 */
	private String invasionStep(Random rng) {
		// get the correct people.
		Person forceOpinion = people.get(rng.nextInt(people.size()));
		ArrayList<Person> neighbors = forceOpinion.neighbors();
		Person changeOpinion = neighbors.get(rng.nextInt(neighbors.size()));
		// see if the forced opinion comes from small group.
		if (forceOpinion.getID() <= smallGroupAmount) {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[0]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[0]--;
			}
			// see if the forced opinion comes from large group.
		} else {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[1]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[1]--;
			}
		}

		changeOpinion.setOpinion(forceOpinion.getOpinion());
//		stateList.add((double) currState[0] / largeGroupAmount + "\t" + (double) currState[1] / smallGroupAmount);
		return forceOpinion.getID() + " -> " + changeOpinion.getID();
	}

	/**
	 * returns the string of the graph which is the list of states.
	 */
	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		for (String s : stateList) {
			str.append(s + "\n");
		}
		return str.toString();
	}

	/**
	 * Returns the current state of the model. The amount of 1 opinions in the large
	 * group followed by the number of 1 opinions in the small group.
	 */
	public String currentState() {
		return currState[0] + "\t" + currState[1];
	}

	/**
	 * Determines whether the model has reached a consensus.
	 */
	public boolean consensus() {
		return (currState[0] == 0 && currState[1] == 0)
				|| (currState[0] == largeGroupAmount && currState[1] == smallGroupAmount);
	}

	/**
	 * Computes the average value of integers in a list.
	 */
	private static double average(ArrayList<Integer> arr) {
		int count = 0;
		for (Integer i : arr) {
			count += i;
		}
		return (double) count / arr.size();
	}

	/**
	 * Computes the sample variance of integers in a list.
	 */
	private static double sampleVariance(ArrayList<Integer> arr) {
		double average = average(arr);
		double variance = 0;
		for (int i = 0; i < arr.size(); i++) {
			variance += Math.pow(arr.get(i) - average, 2) / (arr.size() - 1);
		}
		return variance;
	}

	/**
	 * Creates a hashmap consisting of a string version of each state (i.e. "[5,
	 * 2]").
	 * 
	 * @param largeGroup - amount in large group.
	 * @param smallGroup - amount in small group.
	 */
	private static HashMap<String, Integer> hashMap(int largeGroup, int smallGroup) {

		HashMap<String, Integer> map = new HashMap<>();
		for (int i = 0; i <= Math.max(largeGroup, smallGroup); i++) {
			for (int j = 0; j <= Math.min(largeGroup, smallGroup); j++) {
				map.put(i + "\t" + j, 0);
			}
		}
		return map;
	}

	/*
	 * Converts string into the actual array representation of the state.
	 */
	private static int[] parse(String vect) {
		int[] nums = new int[2];
		int cutoff = vect.indexOf('\t');
		nums[0] = Integer.parseInt(vect.substring(0, cutoff));
		nums[1] = Integer.parseInt(vect.substring(cutoff + 1, vect.length()));
		return nums;
	}

	/**
	 * Creates ordering among the states. Used for sorting to display in a
	 * reasonable fashion.
	 * 
	 * @param s1 - state 1.
	 * @param s2 - state 2.
	 */
	private static int compare(String s1, String s2) {
		int cutoff1 = s1.indexOf('\t');
		int firstEntry1 = Integer.parseInt(s1.substring(0, cutoff1));
		int secondEntry1 = Integer.parseInt(s1.substring(cutoff1 + 1, s1.length()));
		int cutoff2 = s2.indexOf('\t');
		int firstEntry2 = Integer.parseInt(s2.substring(0, cutoff2));
		int secondEntry2 = Integer.parseInt(s2.substring(cutoff2 + 1, s2.length()));

		if (firstEntry1 > firstEntry2) {
			return 1;
		} else if (firstEntry1 < firstEntry2) {
			return -1;
		} else {
			if (secondEntry1 > secondEntry2) {
				return 1;
			} else if (secondEntry1 < secondEntry2) {
				return -1;
			} else {
				return 0;
			}
		}
	}

	/**
	 * Create a new startpoint for the graph by sampling from the distribution of
	 * the proportion of time spent at each state.. Here we define a mult set to be
	 * a set of identical numbers (usually the max of 2147483647). This is due to
	 * java's integer limitation and we use this to do samples beyond that
	 * limitation.
	 * 
	 * @param map   - hashmap consisting of the count for the amount of times each
	 *              state has been reached overall.
	 * @param count - amount of steps for the mult set.
	 * @param step  - total steps per mult set.
	 * @param mult  - amount of mult sets.
	 */
	private int[] startPoint(HashMap<String, Integer> map, int count, int step, int mult) {
		double randVal = new Random().nextDouble();
		double currentInc = 0;
		for (Entry<String, Integer> str : map.entrySet()) {
			// never been visited.
			if (str.getValue() == 0) {
				continue;
			}
			/*
			 * To prevent integer overflow. We are solving for the inequality
			 * currTotal/(step*mult + count) < randVal where currTotal is the cumulative
			 * amount of times at the state for each state checked. We manipulate the
			 * equation be taking the inverse of the left side, separating the numerator,
			 * and spliting the currTotal into the product of square roots. Then we separate
			 * the numerator and denominator to avoid integer overflow.
			 */
			double sOt = (double) step / Math.sqrt(str.getValue());
			double mOt = (double) (mult - 1) / Math.sqrt(str.getValue());
			// if this inequality holds, then it means the amount of steps at the state is
			// so small that we just consider it to be 0.
			if (2147483747.0 / sOt < mOt) {
				continue;
			}
			// converting back into a proportion.
			double cOt = (double) count / str.getValue();
			currentInc += Math.pow(sOt * mOt + cOt, -1);
			if (currentInc > randVal) {
				currState = parse(str.getKey());
				return currState;
			}
		}
		// if no state is chosen. This only happens when randVal is extremely close to
		// 1. We will just restart the process.
		return startPoint(map, count, step, mult);
	}

	/**
	 * Taking multiple random integers in a range without replacement.
	 * 
	 * @param low    - lower bound.
	 * @param high   - upper bound.
	 * @param amount - the amount of numbers to generate
	 */
	private static int[] multipleRandomNumbers(int low, int high, int amount) {
		if (high - low + 1 < amount) {
			amount = high - low + 1;
		}
		Random rng = new Random();
		int[] nums = new int[amount];
		ArrayList<Integer> arr = new ArrayList<>();
		for (int i = low; i <= high; i++) {
			arr.add(i);
		}
		for (int i = 0; i < amount; i++) {
			int rand = rng.nextInt(arr.size());
			nums[i] = arr.get(rand);
			arr.remove(rand);
		}
		return nums;
	}

	/**
	 * We focus on two people and each time the opinion of one of those people is
	 * taken, we move our focus to who gave that opinion. This continues the two
	 * focuses meet together or (coalesce). The thing to observe here is that the
	 * eigenvalue of this transition model is the exact same as the voter model.
	 * 
	 * @param n    - size of the large group.
	 * @param m    - size of the small group.
	 * @param type - 1 = one in each group, 2 = both start in large group, 3 = both
	 *             start in small group, otherwise it is completely random.
	 * @return amount of time to coalesce.
	 */
	public static int coalesceTime(int n, int m, int type) {
		Random rng = new Random();
		int count = 0;
		ArrayList<Person> people = new ArrayList<>();
		for (int i = 1; i <= n + m; i++) {
			people.add(new Person(i));
		}
		// make relations.
		for (int i = 0; i < m; i++) {
			for (int j = m; j < n + m; j++) {
				people.get(i).addNeighbor(people.get(j));
			}
		}
		// make starting points.
		int[] randNums;
		switch (type) {
		case 1:
			randNums = new int[2];
			// one from each.
			randNums[0] = rng.nextInt(m);
			randNums[1] = m + rng.nextInt(n);
			break;
		// both in large group.
		case 2:
			randNums = multipleRandomNumbers(m, n + m - 1, 2);
			break;
		// both in small group.
		case 3:
			if (m > 1) {
				randNums = multipleRandomNumbers(0, m - 1, 2);
			} else {
				randNums = multipleRandomNumbers(0, m + n - 1, 2);
			}
			break;
		// choose randomly if it is not specified.
		default:
			randNums = multipleRandomNumbers(0, m + n - 1, 2);
			break;
		}
		Person p1 = people.get(randNums[0]);
		Person p2 = people.get(randNums[1]);
		// start the simulation.
		while (p1 != p2) {
			count++;
			Person forceOp = people.get(rng.nextInt(people.size()));
			ArrayList<Person> neighbors = forceOp.neighbors();
			Person getOp = neighbors.get(rng.nextInt(neighbors.size()));
			if (getOp == p1) {
				p1 = forceOp;
			} else if (getOp == p2) {
				p2 = forceOp;
			}
		}
		return count;
	}

	/**
	 * Obtain a list for coalesce times.
	 * 
	 * @param n      - size of the large group.
	 * @param m      - size of the small group.
	 * @param type   - 1 = one in each group, 2 = both start in large group, 3 =
	 *               both start in small group, otherwise it is completely random.
	 * @param trials - amount of times to run experiment.
	 */
	public static ArrayList<Integer> coalesceTimeList(int n, int m, int type, int trials) {
		ArrayList<Integer> arr = new ArrayList<>();
		for (int i = 0; i < trials; i++) {
			arr.add(coalesceTime(n, m, type));
		}
		return arr;
	}

	/**
	 * Simulates the quasi-stationary distribution or the proportion of time spent
	 * at each state through the process conditioned on not going to consensus.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param voter - voter model if true. Invasion model if false.
	 * @return
	 */
	public static String quasiStationaryDistribution(int n, int m, int steps, boolean voter) {
		return quasiStationaryDistribution(n, m, steps, 1, voter);
	}

	/**
	 * Simulates the quasi-stationary distribution or the proportion of time spent
	 * at each state through the process conditioned on not going to consensus.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param mult  - amount of times to go "steps" amount of steps (i.e. mult*steps
	 *              total). This is to avoid integer overflow errors.
	 * @param voter - voter model if true. Invasion model if false.
	 * @return
	 */
	public static String quasiStationaryDistribution(int n, int m, int steps, int mult, boolean voter) {
		HashMap<String, Integer> space = quasiStationaryDistributionHash(n, m, steps, mult, voter);
		return convertHashMapToString(n, m, steps, mult, space);
	}

	/**
	 * Given a hashmap, converts into a string representation.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param mult  - amount of times to go "steps" amount of steps (i.e. mult*steps
	 *              total). This is to avoid integer overflow errors.
	 * @param space - the hashmap of all the states.
	 * @return
	 */
	private static String convertHashMapToString(int n, int m, int steps, int mult, HashMap<String, Integer> space) {
		ArrayList<String> arr = new ArrayList<>();
		for (Entry<String, Integer> str : space.entrySet()) {
			arr.add(str.getKey());
		}
		arr.sort((x, y) -> compare(x, y));
		// do not count absorbing or unreachable states.
		arr.remove("0\t0");
		arr.remove(n + "\t" + m);
		arr.remove("0\t" + m);
		arr.remove(n + "\t0");
		StringBuilder s = new StringBuilder();
		for (String str : arr) {
			s.append(str + "\t" + space.get(str) + "\t" + (double) ((double) space.get(str) / steps) / mult + "\n");
		}
		return s.toString();
	}

	/**
	 * Simulates the quasi-stationary distribution or the proportion of time spent
	 * at each state through the process conditioned on not going to consensus.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param mult  - amount of times to go "steps" amount of steps (i.e. mult*steps
	 *              total). This is to avoid integer overflow errors.
	 * @param map   - a hashmap of values.
	 * @return
	 */
	private static String quasiStationaryDistribution(int n, int m, int steps, int mult, HashMap<String, Integer> map) {
		return convertHashMapToString(n, m, steps, mult, map);
	}

	/**
	 * Returns the string representation of the hashmap of the marginal distribution
	 * with the proportion of time spent at each N value.
	 * 
	 * @param n   - number of nodes in big group.
	 * @param map - map of all the amounts of steps for each state.
	 */
	public static String marginalN(int n, HashMap<String, Integer> map) {
		HashMap<Integer, Integer> marg = marginalNWork(n, map);
		StringBuilder s = new StringBuilder();
		int total = total(marg);
		for (Entry<Integer, Integer> num : marg.entrySet()) {
			s.append(num.getKey() + "\t" + num.getValue() + "\t" + (double) num.getValue() / total + "\n");
		}
		return s.toString();
	}

	/**
	 * Returns the hashmap of the marginal distribution with the proportion of time
	 * spent at each N value.
	 * 
	 * @param n   - number of nodes in big group.
	 * @param map - map of all the amounts of steps for each state.
	 */
	private static HashMap<Integer, Integer> marginalNWork(int n, HashMap<String, Integer> map) {
		HashMap<Integer, Integer> marg = new HashMap<>();
		for (int i = 0; i <= n; i++) {
			marg.put(i, 0);
		}
		for (Entry<String, Integer> str : map.entrySet()) {
			String s = str.getKey().substring(0, str.getKey().indexOf('\t'));
			int k = Integer.parseInt(s);
			int amount = marg.get(k) + str.getValue();
			marg.put(k, amount);
		}
		return marg;
	}

	/**
	 * Returns the hashmap of the marginal distribution with the proportion of time
	 * spent at each N value.
	 * 
	 * @param n   - number of nodes in big group.
	 * @param map - map of all the amounts of steps for each state.
	 */
	public static HashMap<Integer, Double> marginalNProp(int n, HashMap<String, Integer> map) {
		HashMap<Integer, Integer> marg = marginalNWork(n, map);
		// to keep track of the proportion of time in each N.
		HashMap<Integer, Double> dub = new HashMap<>();
		// adding all elements in.
		for (int i = 0; i <= n; i++) {
			dub.put(i, 0.0);
		}
		int total = total(marg);
		// proportions for the marginal.
		for (Entry<Integer, Integer> num : marg.entrySet()) {
			dub.put(num.getKey(), (double) num.getValue() / total);
		}
		return dub;
	}

	/**
	 * Returns the hashmap of the marginal distribution with the proportion of time
	 * spent at each N value.
	 * 
	 * @param n    - number of nodes in big group.
	 * @param marg - map of all the amounts of steps for each N.
	 */
	private static HashMap<Integer, Double> marginalNPropInputMarginal(int n, HashMap<Integer, Integer> marg) {
		// to keep track of the proportion of time in each N.
		HashMap<Integer, Double> dub = new HashMap<>();
		// adding all elements in.
		for (int i = 0; i <= n; i++) {
			dub.put(i, 0.0);
		}
		int total = total(marg);
		// proportions for the marginal.
		for (Entry<Integer, Integer> num : marg.entrySet()) {
			dub.put(num.getKey(), (double) num.getValue() / total);
		}
		return dub;
	}

	/**
	 * Iterates through each i = 0,...,n with the given probability of being at that
	 * state from the weights variable. Then computes the binomial distribution at
	 * each of those points which results in a mixture.
	 * 
	 * @param n   - number of nodes in the big group.
	 * @param m   - number of nodes in the small group.
	 * @param map - the hashmap of the marginal distribution of the amount in the
	 *            big group.
	 */
	public static HashMap<Integer, Double> mixtureN(int n, int m, HashMap<Integer, Double> weights) {
		// keeps track of mixture probabilities.
		HashMap<Integer, Double> prob = new HashMap<>();
		for (int i = 0; i <= m; i++) {
			prob.put(i, 0.0);
		}
		// loop through each possibility with its weight.
		for (Entry<Integer, Double> num : weights.entrySet()) {
			BinomialDistribution b = new BinomialDistribution(m, (double) num.getKey() / n);
			// loop through each binomial outcome.
			for (int i = 0; i <= m; i++) {
				// add the probability * weight to what is already there.
				prob.put(i, prob.get(i) + num.getValue() * b.probability(i));
			}
		}
		return prob;
	}

	/**
	 * Iterates through each i = 0,...,n with the given probability of being at that
	 * state from the weights variable. Then computes the binomial distribution at
	 * each of those points which results in a mixture.
	 * 
	 * @param n   - number of nodes in the big group.
	 * @param m   - number of nodes in the small group.
	 * @param map - the hashmap of the marginal distribution of the amount in the
	 *            big group.
	 */
	public static String mixtureNToString(int n, int m, HashMap<Integer, Double> weights) {
		HashMap<Integer, Double> prob = mixtureN(n, m, weights);
		StringBuilder s = new StringBuilder("\n");
		for (Entry<Integer, Double> num : prob.entrySet()) {
			s.append(num.getKey() + "\t" + num.getValue() + "\n");
		}
		return s.toString();
	}

	/**
	 * Returns the string representation of the hashmap of the marginal distribution
	 * with the proportion of time spent at each M value.
	 * 
	 * @param m   - number of nodes in small group.
	 * @param map - map of all the amounts of steps for each state.
	 */
	public static String marginalM(int m, HashMap<String, Integer> map) {
		HashMap<Integer, Integer> marg = new HashMap<>();
		// for each M value, make a map computing the total.
		for (int i = 0; i <= m; i++) {
			marg.put(i, 0);
		}
		// get the key to know where to increment.
		for (Entry<String, Integer> str : map.entrySet()) {
			String s = str.getKey().substring(str.getKey().indexOf('\t') + 1, str.getKey().length());
			int k = Integer.parseInt(s);
			int amount = marg.get(k) + str.getValue();
			marg.put(k, amount);
		}
		// take the total to proportion everything out.
		int total = total(marg);
		StringBuilder s = new StringBuilder();
		for (Entry<Integer, Integer> num : marg.entrySet()) {
			s.append(num.getKey() + "\t" + num.getValue() + "\t" + (double) num.getValue() / total + "\n");
		}
		return s.toString();
	}

	/**
	 * Given each M value, finds the proportion of time spent in N for all M. Gives
	 * a list for each M.
	 * 
	 * @param n   - number of nodes in the big group.
	 * @param m   - number of nodes in the small group.
	 * @param map - the hashmap of the marginal distribution of the amount in the
	 *            big group.
	 */
	public static String conditionalsOfAllN(int n, int m, HashMap<String, Integer> map) {
		HashMap<Integer, HashMap<Integer, Integer>> cond = new HashMap<>();
		for (int i = 0; i <= n; i++) {
			cond.put(i, new HashMap<Integer, Integer>());
		}
		for (Entry<String, Integer> str : map.entrySet()) {
			String parseN = str.getKey().substring(0, str.getKey().indexOf('\t'));
			int nVal = Integer.parseInt(parseN);
			String parseM = str.getKey().substring(str.getKey().indexOf('\t') + 1, str.getKey().length());
			int mVal = Integer.parseInt(parseM);
			cond.get(nVal).put(mVal, str.getValue());
		}
		StringBuilder s = new StringBuilder();
		for (Entry<Integer, HashMap<Integer, Integer>> vals : cond.entrySet()) {
			int total = total(vals.getValue());
			if (total == 0) {
				continue;
			}
			s.append("j=" + vals.getKey() + "\n");
			for (Entry<Integer, Integer> ints : vals.getValue().entrySet()) {
				if ((vals.getKey() == 0 || vals.getKey() == n) && (ints.getKey() == 0 || ints.getKey() == m)) {
					continue;
				}
				s.append(ints.getKey() + "\t" + ints.getValue() + "\t" + (double) ints.getValue() / total + "\n");
			}
			s.append("\n");
		}
		return s.toString();
	}

	/**
	 * Given each N value, finds the proportion of time spent in M for all M. Gives
	 * a list for each N.
	 * 
	 * @param n   - number of nodes in the big group.
	 * @param m   - number of nodes in the small group.
	 * @param map - the hashmap of the marginal distribution of the amount in the
	 *            big group.
	 */
	public static String conditionalsOfAllM(int n, int m, HashMap<String, Integer> map) {
		HashMap<Integer, HashMap<Integer, Integer>> cond = new HashMap<>();
		for (int i = 0; i <= m; i++) {
			cond.put(i, new HashMap<Integer, Integer>());
		}
		for (Entry<String, Integer> str : map.entrySet()) {
			String parseN = str.getKey().substring(0, str.getKey().indexOf('\t'));
			int nVal = Integer.parseInt(parseN);
			String parseM = str.getKey().substring(str.getKey().indexOf('\t') + 1, str.getKey().length());
			int mVal = Integer.parseInt(parseM);
			// setting each m value to have its own list.
			cond.get(mVal).put(nVal, str.getValue());
		}
		StringBuilder s = new StringBuilder();
		for (Entry<Integer, HashMap<Integer, Integer>> vals : cond.entrySet()) {
			int total = total(vals.getValue());
			if (total == 0) {
				continue;
			}
			s.append("i=" + vals.getKey() + "\n");
			for (Entry<Integer, Integer> ints : vals.getValue().entrySet()) {
				// take out absorbing states and unreachable states instead of having them show
				// up as 0.
				if ((vals.getKey() == 0 || vals.getKey() == m) && (ints.getKey() == 0 || ints.getKey() == n)) {
					continue;
				}
				s.append(ints.getKey() + "\t" + ints.getValue() + "\t" + (double) ints.getValue() / total + "\n");
			}
			s.append("\n");
		}
		return s.toString();
	}

	/**
	 * Simulates the quasi-stationary distribution.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param voter - voter model if true. Invasion model if false.
	 */
	private static HashMap<String, Integer> quasiStationaryDistributionHash(int n, int m, int steps, boolean voter) {
		return quasiStationaryDistributionHash(n, m, steps, 1, voter);
	}

	/**
	 * Determines whether there should be a progress update for how long the code
	 * has been running. Usually in 10% increments.
	 * 
	 * @param currStep - current step count.
	 * @param currMult - current mult count.
	 * @param steps    - total steps.
	 * @param mult     - total amount of sets of "steps".
	 * @return true if update, false otherwise.
	 */
	private static boolean update(int currStep, int currMult, int steps, int mult) {
		// mult >= 10 case.
		if (mult >= 10) {
			if ((currMult % (mult / 10) == 0 || currMult == mult) && currStep == steps) {
				return true;
			}
			return false;
		}
		// mult < 10 case.
		int times = 10 / mult;
		if ((steps > times - 1 && currStep % (steps / times) == 0) || currStep == steps) {
			return true;
		}
		return false;
	}

	/**
	 * Simulates the quasi-stationary distribution.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param mult  - total amount of sets of "steps".
	 * @param voter - voter model if true. Invasion model if false.
	 */
	private static HashMap<String, Integer> quasiStationaryDistributionHash(int n, int m, int steps, int mult,
			boolean voter) {
		HashMap<String, Integer> space = hashMap(n, m);

		// counters.
		int curr = 0;
		int multCounter = 1;
		// keep status for how many times we reach consensus.
		int graphMakeCounter = -1;
		int graphTotal = 0;
		BipartiteGraph graph = new BipartiteGraph(n, m, voter);
		Random rng = new Random();
		// choose random from 1 to n-1 then random from 1 to m to guarantee no consensus
		// to start.
		int rand1 = rng.nextInt(n - 1) + 1;
		int rand2 = rng.nextInt(m) + 1;
		graph.set(rand1, rand2);

		// keep reference for how long how code will run.
		System.out.println("Start");

		long startTime = System.currentTimeMillis();
		while (curr < steps) {
			// reset the amount of steps on the graph so we don't have to make a new graph.
			graph.reset();
			// if it is not the first step.
			if (curr != 0 || multCounter != 1) {
				int[] state = graph.startPoint(space, curr, steps, multCounter);
				if (state != null) {
					graph.set(state[0], state[1]);
				}
			}
			// re-randomize everything if we initialize a state that is unreachable or
			// absorbing.
			if (graph.consensus() || graph.currentState().contains("0\t" + m)
					|| graph.currentState().contains(n + "\t0")) {
				continue;
			}
			graphMakeCounter++;
			// keep following these steps until consensus or number of steps.
			while (!graph.consensus() && curr < steps) {
				curr++;
				// update periodically so we know the progress of the code.
				if (update(curr, multCounter, steps, mult)) {
					graphTotal += graphMakeCounter;
					System.out.println(multCounter + "\t" + curr + "   (" + n + "," + m
							+ ")\tAmount of Consensus since last update: " + graphMakeCounter
							+ "\tConsensus amount total: " + graphTotal + "\tCurrent State: " + graph.currentState());
					graphMakeCounter = 0;
				}
				// reset the step counter.
				if (curr == steps && multCounter < mult) {
					curr = 0;
					multCounter++;
				}
				// change the state and keep the data.
				String currentState = graph.currentState();
				int amount = space.get(currentState) + 1;
				space.put(currentState, amount);
				graph.step();
			}
		}
		long endTime = System.currentTimeMillis();
		System.out.println(endTime - startTime);
		return space;
	}

	/**
	 * Simulates the quasi-stationary distribution in table format.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param voter - voter model if true. Invasion model if false.
	 */
	public static String quasiStationaryDistributionTable(int n, int m, int steps, boolean voter) {
		return quasiStationaryDistributionTable(n, m, steps, 1, voter);
	}

	/**
	 * Simulates the quasi-stationary distribution in table format.
	 * 
	 * @param n     - size of the large group.
	 * @param m     - size of the small group.
	 * @param steps - amount of steps to run simulation for.
	 * @param mult  - total amount of sets of "steps".
	 * @param voter - voter model if true. Invasion model if false.
	 */
	public static String quasiStationaryDistributionTable(int n, int m, int steps, int mult, boolean voter) {
		StringBuilder s = new StringBuilder();
		HashMap<String, Integer> space = hashMap(n, m);
		// keep track of the N then M value for the nested hashmap.
		HashMap<Integer, HashMap<Integer, Double>> map = new HashMap<>();
		for (int i = 0; i <= n; i++) {
			HashMap<Integer, Double> dub = new HashMap<>();
			map.put(i, dub);
			for (int j = 0; j <= m; j++) {
				dub.put(j, 0.0);
			}
		}
		// counters.
		int curr = 0;
		int multCounter = 1;
		// keep status for how many times we reach consensus.
		int graphMakeCounter = -1;
		int graphTotal = 0;
		BipartiteGraph graph = new BipartiteGraph(n, m, voter);
		Random rng = new Random();
		// choose random from 1 to n-1 then random from 1 to m to guarantee no consensus
		// to start.
		int rand1 = rng.nextInt(n - 1) + 1;
		int rand2 = rng.nextInt(m) + 1;
		graph.set(rand1, rand2);

		// keep reference for how long how code will run.
		System.out.println("Start");

		long startTime = System.currentTimeMillis();
		while (curr < steps) {
			// reset the amount of steps on the graph so we don't have to make a new graph.
			graph.reset();
			// if it is not the first step.
			if (curr != 0 || multCounter != 1) {
				int[] state = graph.startPoint(space, curr, steps, multCounter);
				if (state != null) {
					graph.set(state[0], state[1]);
				}
			}
			// re-randomize everything if we initialize a state that is unreachable or
			// absorbing.
			if (graph.consensus() || graph.currentState().contains("0\t" + m)
					|| graph.currentState().contains(n + "\t0")) {
				continue;
			}
			graphMakeCounter++;
			// keep following these steps until consensus or number of steps.
			while (!graph.consensus() && curr < steps) {
				curr++;
				// update periodically so we know the progress of the code.
				if (update(curr, multCounter, steps, mult)) {
					graphTotal += graphMakeCounter;
					System.out.println(multCounter + "\t" + curr + "   (" + n + "," + m
							+ ")\tAmount of Consensus since last update: " + graphMakeCounter
							+ "\tConsensus amount total: " + graphTotal + "\tCurrent State: " + graph.currentState());
					graphMakeCounter = 0;
				}
				// reset the step counter.
				if (curr == steps && multCounter < mult) {
					curr = 0;
					multCounter++;
				}
				// change the state and keep the data.
				String currentState = graph.currentState();
				int amount = space.get(currentState) + 1;
				space.put(currentState, amount);
				graph.step();
				// put into list
				String parsJ = currentState.substring(0, currentState.indexOf('\t'));
				int amountJ = Integer.parseInt(parsJ);
				String parsI = currentState.substring(currentState.indexOf('\t') + 1, currentState.length());
				int amountI = Integer.parseInt(parsI);
				Double d = map.get(amountJ).get(amountI);
				d += (double) ((double) 1 / steps) / mult;
				map.get(amountJ).put(amountI, d);
			}
		}
		long endTime = System.currentTimeMillis();
		System.out.println(endTime - startTime);
		// output the data.
		for (Entry<Integer, HashMap<Integer, Double>> ent : map.entrySet()) {
			for (Entry<Integer, Double> nums : ent.getValue().entrySet()) {
				s.append(nums.getValue() + "\t");
			}
			s.append("\n");
		}
		return s.toString();
	}

	/**
	 * Average time to coalesce for a range of 3,...,n and 1,...,m values.
	 * 
	 * @param n - upper bound for number of people in the large group.
	 * @param m - upper bound for number of people in the small group.
	 */
	public static String coalesceTimeData(int n, int m) {
		int trials = 1000;
		StringBuilder output = new StringBuilder("(n,m) ");
		// nested array list to represent m then n values.
		ArrayList<ArrayList<Double>> arr = new ArrayList<>();
		// input the m value columns.
		for (int j = 1; j <= m; j++) {
			output.append(j + "\t");
			ArrayList<Double> averages = new ArrayList<>();
			arr.add(averages);
			// get data for each corresponding n.
			for (int i = 3; i <= n; i++) {
				averages.add(average(coalesceTimeList(i, j, 0, trials)));
				// update the progress.
				if ((n > 3 && i % (n / 4) == 0) || i == n) {
					System.out.println(j + "\t" + i);
				}
			}
		}
		System.out.println("printing");
		output.append("\n");
		// output the data.
		for (int i = 0; i < arr.get(0).size(); i++) {
			output.append((i + 3) + "\t");
			for (int j = 0; j < arr.size(); j++) {
				output.append(arr.get(j).get(i) + "\t");
			}
			output.append("\n");
		}
		System.out.println("done");
		return output.toString();
	}

	/**
	 * Average time to coalesce with a range of n values from 3,...,n and m = 1
	 * 
	 * @param n      - upper bound for number of people in the large group.
	 * @param trials - number of trials for each.
	 */
	public static String coalesceHistograms(int n, int trials) {
		StringBuilder output = new StringBuilder();
		ArrayList<ArrayList<Integer>> arr = new ArrayList<>();
		for (int i = 3; i <= n; i++) {
			output.append(i + " ");
			arr.add(coalesceTimeList(i, 1, 0, trials));
			// give status update.
			if ((n > 3 && i % (n / 4) == 0) || i == n) {
				System.out.println(i);
			}
		}
		System.out.println("printing");
		output.append("\n");
		for (int i = 0; i < trials; i++) {
			for (int j = 0; j < arr.size(); j++) {
				output.append(arr.get(j).get(i) + " ");
			}
			output.append("\n");
		}
		System.out.println("done");
		return output.toString();
	}

	/**
	 * Computes the total integer value amount all keys in the hashmap.
	 * 
	 * @param map - Hashmap of integer keys and integer values.
	 * @return the total.
	 */
	private static int total(HashMap<Integer, Integer> map) {
		int total = 0;
		for (Entry<Integer, Integer> str : map.entrySet()) {
			total += str.getValue();
		}
		return total;
	}

	/**
	 * Outputs all the useful data such as the proportion of time spent at each
	 * state, the marginals, the and the conditionals.
	 * 
	 * @param n     - size of big group
	 * @param m     - size of small group
	 * @param steps - number of steps
	 * @param voter - voter model or invasion model.
	 */
	public static String allQSDData(int n, int m, int steps, boolean voter) {
		return allQSDData(n, m, steps, 1, voter);
	}

	/**
	 * Outputs all the useful data such as the proportion of time spent at each
	 * state, the marginals, the and the conditionals.
	 * 
	 * @param n     - size of big group
	 * @param m     - size of small group
	 * @param steps - number of steps
	 * @param mult  - multiplicity of amount of steps (to go beyond the limit)
	 * @param voter - voter model or invasion model.
	 */
	public static String allQSDData(int n, int m, int steps, int mult, boolean voter) {
		StringBuilder s = new StringBuilder("");
		// generate qsd data.
		HashMap<String, Integer> map = quasiStationaryDistributionHash(n, m, steps, mult, voter);
		// generate marginal data
		HashMap<Integer, Integer> marg = marginalNWork(n, map);
		// generate proportions for each N.
		HashMap<Integer, Double> prop = marginalNPropInputMarginal(n, marg);
		s.append("Quasi-Stationary Distribution\n" + quasiStationaryDistribution(n, m, steps, mult, map));
		s.append("\nMarginal Distribution of N\n" + marginalN(n, map));
		s.append("\nMarginal Distribution of M\n" + marginalM(m, map));
//		s.append("\nBinomial Mixture M" + mixtureNToString(n, m, prop));
		s.append("\nConditional Distribution M|N\n" + conditionalsOfAllN(n, m, map));
		s.append("\nConditional Distribution N|M\n" + conditionalsOfAllM(n, m, map));
		return s.toString();
	}

	/**
	 * Average amount of time for the model to come to consensus given a set start
	 * point.
	 * 
	 * @param amount - amount of times to run the experiment.
	 * @param n      - number of people in large group.
	 * @param m      - number of people in small group.
	 * @param state  - starting state.
	 * @param voter  - voter model if true. Invasion model if false.
	 */
	public static double timeToAbsorbtion(int amount, int n, int m, String state, boolean voter) {
		int[] pars = parse(state);
		int count = 0;
		double avg = 0.0;
		// running the experiment.
		while (count < amount) {
			BipartiteGraph graph = new BipartiteGraph(n, m, voter);
			graph.set(pars[0], pars[1]);
			graph.simulation();
			avg += (double) graph.time() / amount;
			count++;
		}
		return avg;
	}

	public static void main(String[] args) {
		try {
			PrintWriter out = new PrintWriter("src/voterInvasionModels/Invasion.txt");
			int n = 100;
			int m = 10;
			int steps = Integer.MAX_VALUE;
			int mult = 5;
			out.println(quasiStationaryDistributionTable(n, m, steps, mult, true));
//			out.println(allQSDData(n, m, 100000, 10, true));
//			BipartiteGraph g = new BipartiteGraph(200, 1, false);
//			g.simulation();
//			out.println(g.toString());
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
