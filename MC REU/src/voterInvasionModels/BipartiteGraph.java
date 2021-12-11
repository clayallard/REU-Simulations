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
 * @author Clayton Allard
 *
 */
public class BipartiteGraph {

	private ArrayList<Person> people = new ArrayList<>();
//	private HashSet<String> connections = new HashSet<>();
	private boolean voterModel;
	private int time;
	private int timeToAbsorbing;
	private int smallGroupAmount;
	private int largeGroupAmount;
	private int[] currState = new int[2];
	private ArrayList<String> stateList = new ArrayList<>();

	public void reset() {
		time = 0;
		timeToAbsorbing = 0;
		stateList = new ArrayList<>();
	}

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

	private void addRelation(Person p1, Person p2) {
		if (p1.equals(p2)) {
			return;
		}
		p1.addNeighbor(p2);
		String s = connectionString(p1, p2);
//		if (!connections.contains(s)) {
//			connections.add(s);
//		}
	}

	private void set(int k, int h) {
		currState[0] = k;
		currState[1] = h;
		for (int i = 0; i < smallGroupAmount; i++) {
			if (i < h) {
				people.get(i).setOpinion(1);
			} else {
				people.get(i).setOpinion(0);
			}
		}
		for (int i = smallGroupAmount; i < people.size(); i++) {
			if (i < k + smallGroupAmount) {
				people.get(i).setOpinion(1);
			} else {
				people.get(i).setOpinion(0);
			}
		}
		stateList.add((double) currState[0] / largeGroupAmount + "\t" + (double) currState[1] / smallGroupAmount);
	}

	public void addRelations(ArrayList<Person> p1, ArrayList<Person> p2) throws IllegalArgumentException {
		if (p1.size() != p2.size()) {
			throw new IllegalArgumentException("Must have equal sized lists");
		}
		// don't want to change anything until it is safe to.
		for (int i = 0; i < p1.size(); i++) {
			if (p1.get(i).equals(p2.get(i))) {
				throw new IllegalArgumentException("Cannot have self relation");
			}
		}
		for (int i = 0; i < p1.size(); i++) {
			addRelation(p1.get(i), p2.get(i));
		}
	}

	public void addRelations(ArrayList<Person> peep) {
		for (int i = 0; i < peep.size() - 1; i++) {
			for (int j = i + 1; j < peep.size(); j++) {
				addRelation(peep.get(i), peep.get(j));
			}
		}
	}

	public ArrayList<Person> allNeighbors(Person p) {
		return p.neighbors();
	}

//	public String connectionString() {
//		String str = "";
//		ArrayList<String> arr = new ArrayList<>(connections);
//		arr.sort(null);
//		for (String s : arr) {
//			str += s + "\n";
//		}
//		return str;
//	}

	private String connectionString(Person p1, Person p2) {

		if (p1.getID() < p2.getID()) {
			return p1.getID() + " -- " + p2.getID();
		} else {
			return p2.getID() + " -- " + p1.getID();
		}
	}

	public void removeRelation(Person p1, Person p2) {
		if (p1.removeNeighbor(p2) && p2.removeNeighbor(p1)) {
//			connections.remove(connectionString(p1, p2));
			if (p1.neighbors().isEmpty()) {
				people.remove(p1);
			}
			if (p2.neighbors().isEmpty()) {
				people.remove(p2);
			}
		}
	}

	public void removeAllRelations(Person p) {
		for (Person per : p.neighbors()) {
			removeRelation(p, per);
		}
	}

	public void removeRelations(ArrayList<Person> p1, ArrayList<Person> p2) throws IllegalArgumentException {
		if (p1.size() != p2.size()) {
			throw new IllegalArgumentException("Must have equal sized lists");
		}
		for (int i = 0; i < p1.size(); i++) {
			removeRelation(p1.get(i), p2.get(i));
		}
	}

	public String step() {
		Random rng = new Random();
		return step(rng);
	}

	public String step(Random rng) {
		time++;
		if (voterModel) {
			return voterStep(rng);
		} else {
			return invasionStep(rng);
		}
	}

	public void simulation() {
		simulation(new Random());
	}

	public void simulation(Random rng) {
		while (!consensus())
			step(rng);
	}

	public int time() {
		return time;
	}

	public int timeToConsensus() {
		return timeToAbsorbing;
	}

	private String voterStep(Random rng) {
		Person changeOpinion = people.get(rng.nextInt(people.size()));
		ArrayList<Person> neighbors = changeOpinion.neighbors();
		Person forceOpinion = neighbors.get(rng.nextInt(neighbors.size()));
		if (forceOpinion.getID() <= smallGroupAmount) {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[0]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[0]--;
			}
		} else {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[1]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[1]--;
			}
		}

		changeOpinion.setOpinion(forceOpinion.getOpinion());
//		stateList.add((double) currState[0] / largeGroupAmount + "\t" + (double) currState[1] / smallGroupAmount);
//		if (currentAmount() != 0 && currentAmount() != people.size()) {
//			timeToAbsorbing = time;
//		}
//		currentAmount.add(currentAmount());
		return changeOpinion.getID() + " <- " + forceOpinion.getID();
	}

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
		} else {
			if (forceOpinion.getOpinion() == 1 && changeOpinion.getOpinion() == 0) {
				currState[1]++;
			} else if (forceOpinion.getOpinion() == 0 && changeOpinion.getOpinion() == 1) {
				currState[1]--;
			}
		}

		changeOpinion.setOpinion(forceOpinion.getOpinion());
//		stateList.add((double) currState[0] / largeGroupAmount + "\t" + (double) currState[1] / smallGroupAmount);
//		if (currentAmount() != 0 && currentAmount() != people.size()) {
//			timeToAbsorbing = time;
//		}
//		currentAmount.add(currentAmount());
		return forceOpinion.getID() + " -> " + changeOpinion.getID();
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		for (String s : stateList) {
			str.append(s + "\n");
		}
		return str.toString();
	}

	public String currentState() {
		return currState[0] + "\t" + currState[1];
	}

	public boolean consensus() {
		return (currState[0] == 0 && currState[1] == 0)
				|| (currState[0] == largeGroupAmount && currState[1] == smallGroupAmount);
	}

	private static double average(ArrayList<Integer> arr) {
		int count = 0;
		for (Integer i : arr) {
			count += i;
		}
		return (double) count / arr.size();
	}

	private static double sampleVariance(ArrayList<Integer> arr) {
		double average = average(arr);
		double variance = 0;
		for (int i = 0; i < arr.size(); i++) {
			variance += Math.pow(arr.get(i) - average, 2) / (arr.size() - 1);
		}
		return variance;
	}

	private static HashMap<String, Integer> hashMap(int largeGroup, int smallGroup) {

		HashMap<String, Integer> map = new HashMap<>();
		for (int i = 0; i <= Math.max(largeGroup, smallGroup); i++) {
			for (int j = 0; j <= Math.min(largeGroup, smallGroup); j++) {
				map.put(i + "\t" + j, 0);
			}
		}
		return map;
	}

	private static int[] parse(String vect) {
		int[] nums = new int[2];
		int cutoff = vect.indexOf('\t');
		nums[0] = Integer.parseInt(vect.substring(0, cutoff));
		nums[1] = Integer.parseInt(vect.substring(cutoff + 1, vect.length()));
		return nums;
	}

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

	private int[] startPoint(HashMap<String, Integer> map, int count, int step, int mult) {
		double randVal = new Random().nextDouble();
		double currentInc = 0;
		for (Entry<String, Integer> str : map.entrySet()) {
			if (str.getValue() == 0) {
				continue;
			}
			double sOt = (double) step / Math.sqrt(str.getValue());
			double mOt = (double) (mult - 1) / Math.sqrt(str.getValue());
			if (2147483747.0 / sOt < mOt) {
				continue;
			}
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
		// start the simultion.
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

	public static ArrayList<Integer> coalesceTimeList(int n, int m, int type, int trials) {
		ArrayList<Integer> arr = new ArrayList<>();
		for (int i = 0; i < trials; i++) {
			arr.add(coalesceTime(n, m, type));
		}
		return arr;
	}

	public static String quasiStationaryDistribution(int n, int m, int steps, boolean voter) {
		return quasiStationaryDistribution(n, m, steps, 1, voter);
	}

	public static String quasiStationaryDistribution(int n, int m, int steps, int mult, boolean voter) {
		HashMap<String, Integer> space = quasiStationaryDistributionHash(n, m, steps, mult, voter);
		return convertHashMapToString(n, m, steps, mult, space);
	}

	private static String convertHashMapToString(int n, int m, int steps, int mult, HashMap<String, Integer> space) {
		ArrayList<String> arr = new ArrayList<>();
		for (Entry<String, Integer> str : space.entrySet()) {
			arr.add(str.getKey());
		}
		arr.sort((x, y) -> compare(x, y));
		arr.remove("0\t0");
		arr.remove(n + "\t" + m);
		arr.remove("0\t" + m);
		arr.remove(n + "\t0");
		String s = "";
		for (String str : arr) {
			s += str + "\t" + space.get(str) + "\t" + (double) ((double) space.get(str) / steps) / mult + "\n";
		}
		return s;
	}

	public static String quasiStationaryDistribution(int n, int m, int steps, int mult, HashMap<String, Integer> map) {
		return convertHashMapToString(n, m, steps, mult, map);
	}

	public static String marginalN(int n, HashMap<String, Integer> map) {
		HashMap<Integer, Integer> marg = marginalNWork(n, map);
		String s = "";
		int total = total(marg);
		for (Entry<Integer, Integer> num : marg.entrySet()) {
			s += num.getKey() + "\t" + num.getValue() + "\t" + (double) num.getValue() / total + "\n";
		}
		return s;
	}

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
	 * @param map - the hashmap of the marginal distribution of the amount in the
	 *            big group.
	 */
	public static String mixtureNToString(int n, int m, HashMap<Integer, Double> weights) {
		HashMap<Integer, Double> prob = mixtureN(n, m, weights);
		// string to add to
		StringBuilder s = new StringBuilder("\n");
		for (Entry<Integer, Double> num : prob.entrySet()) {
			s.append(num.getKey() + "\t" + num.getValue() + "\n");
		}
		return s.toString();
	}

	public static String marginalM(int m, HashMap<String, Integer> map) {
		HashMap<Integer, Integer> marg = new HashMap<>();
		for (int i = 0; i <= m; i++) {
			marg.put(i, 0);
		}
		for (Entry<String, Integer> str : map.entrySet()) {
			String s = str.getKey().substring(str.getKey().indexOf('\t') + 1, str.getKey().length());
			int k = Integer.parseInt(s);
			int amount = marg.get(k) + str.getValue();
			marg.put(k, amount);
		}
		int total = total(marg);
		String s = "";
		for (Entry<Integer, Integer> num : marg.entrySet()) {
			s += num.getKey() + "\t" + num.getValue() + "\t" + (double) num.getValue() / total + "\n";
		}
		return s;
	}

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
		String s = "";
		for (Entry<Integer, HashMap<Integer, Integer>> vals : cond.entrySet()) {
			int total = total(vals.getValue());
			if (total == 0) {
				continue;
			}
			s += "j=" + vals.getKey() + "\n";
			for (Entry<Integer, Integer> ints : vals.getValue().entrySet()) {
				if ((vals.getKey() == 0 || vals.getKey() == n) && (ints.getKey() == 0 || ints.getKey() == m)) {
					continue;
				}
				s += ints.getKey() + "\t" + ints.getValue() + "\t" + (double) ints.getValue() / total + "\n";
			}
			s += "\n";
		}
		return s;
	}

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
		String s = "";
		for (Entry<Integer, HashMap<Integer, Integer>> vals : cond.entrySet()) {
			int total = total(vals.getValue());
			if (total == 0) {
				continue;
			}
			s += "i=" + vals.getKey() + "\n";
			for (Entry<Integer, Integer> ints : vals.getValue().entrySet()) {
				// take out absorbing states and unreachable states instead of having them show
				// up as 0.
				if ((vals.getKey() == 0 || vals.getKey() == m) && (ints.getKey() == 0 || ints.getKey() == n)) {
					continue;
				}
				s += ints.getKey() + "\t" + ints.getValue() + "\t" + (double) ints.getValue() / total + "\n";
			}
			s += "\n";
		}
		return s;
	}

	private static HashMap<String, Integer> quasiStationaryDistributionHash(int n, int m, int steps, boolean voter) {
		return quasiStationaryDistributionHash(n, m, steps, 1, voter);
	}

	private static boolean update(int currStep, int currMult, int steps, int mult) {
		if (mult >= 10) {
			if ((currMult % (mult / 10) == 0 || currMult == mult) && currStep == steps) {
				return true;
			}
			return false;
		}
		int times = 10 / mult;
		if ((steps > times - 1 && currStep % (steps / times) == 0) || currStep == steps) {
			return true;
		}
		return false;
	}

	private static HashMap<String, Integer> quasiStationaryDistributionHash(int n, int m, int steps, int mult,
			boolean voter) {
		HashMap<String, Integer> space = hashMap(n, m);

		int curr = 0;
		int multCounter = 1;
		int graphMakeCounter = -1;
		int graphTotal = 0;
		BipartiteGraph graph = new BipartiteGraph(n, m, voter);

		System.out.println("Start");

		long startTime = System.currentTimeMillis();
		while (curr < steps) {
			graph.reset();
			if (curr != 0) {
				int[] state = graph.startPoint(space, curr, steps, multCounter);
				if (state != null) {
					graph.set(state[0], state[1]);
				}
			}
			if (graph.consensus() || graph.currentState().contains("0\t" + m)
					|| graph.currentState().contains(n + "\t0")) {
				continue;
			}
			graphMakeCounter++;
			while (!graph.consensus() && curr < steps) {
				curr++;
				if (update(curr, multCounter, steps, mult)) {
					graphTotal += graphMakeCounter;
					System.out.println(multCounter + "\t" + curr + "   (" + n + "," + m
							+ ")\tAmount of Consensus since last update: " + graphMakeCounter
							+ "\tConsensus amount total: " + graphTotal + "\tCurrent State: " + graph.currentState());
					graphMakeCounter = 0;
				}
				if (curr == steps && multCounter < mult) {
					curr = 0;
					multCounter++;
				}
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

	public static String quasiStationaryDistributionTable(int n, int m, int steps, boolean voter) {
		return quasiStationaryDistributionTable(n, m, steps, 1, voter);
	}

	public static String quasiStationaryDistributionTable(int n, int m, int steps, int mult, boolean voter) {
		String s = "";
		HashMap<String, Integer> space = hashMap(n, m);
		HashMap<Integer, HashMap<Integer, Double>> map = new HashMap<>();
		for (int i = 0; i <= n; i++) {
			HashMap<Integer, Double> dub = new HashMap<>();
			map.put(i, dub);
			for (int j = 0; j <= m; j++) {
				dub.put(j, 0.0);
			}
		}
		int curr = 0;
		int multCounter = 1;
		int graphMakeCounter = 0;
		int graphTotal = 0;
		BipartiteGraph graph = new BipartiteGraph(n, m, voter);
		while (curr < steps) {
			graph.reset();
			if (curr != 0) {
				int[] state = graph.startPoint(space, curr, steps, multCounter);
				if (state != null) {
					graph.set(state[0], state[1]);
				}
			}
			if (graph.consensus() || graph.currentState().contains("0\t" + m)
					|| graph.currentState().contains(n + "\t0")) {
				continue;
			}
			graphMakeCounter++;
			while (!graph.consensus() && curr < steps) {
				curr++;
				if (update(curr, multCounter, steps, mult)) {
					graphTotal += graphMakeCounter;
					System.out.println(multCounter + "\t" + curr + "   (" + n + "," + m
							+ ")\tAmount of Consensus since last update: " + graphMakeCounter
							+ "\tConsensus amount total: " + graphTotal);
					graphMakeCounter = 0;
				}
				if (curr == steps && multCounter < mult) {
					curr = 0;
					multCounter++;
				}
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
		for (Entry<Integer, HashMap<Integer, Double>> ent : map.entrySet()) {
			for (Entry<Integer, Double> nums : ent.getValue().entrySet()) {
				s += nums.getValue() + "\t";
			}
			s += "\n";
		}
		return s;
	}

	public static String coalesceTimeData(int n, int m) {
		int trials = 1000;
		String output = "(n,m) ";
		ArrayList<ArrayList<Double>> arr = new ArrayList<>();
		for (int j = 1; j <= m; j++) {
			output += j + " ";
			ArrayList<Double> averages = new ArrayList<>();
			arr.add(averages);
			for (int i = 3; i <= n; i++) {
				averages.add(average(coalesceTimeList(i, j, 0, trials)));
				if ((n > 3 && i % (n / 4) == 0) || i == n) {
					System.out.println(j + "  " + i);
				}
			}
		}
		System.out.println("printing");
		output += "\n";
		for (int i = 0; i < arr.get(0).size(); i++) {
			output += (i + 3) + " ";
			for (int j = 0; j < arr.size(); j++) {
				output += arr.get(j).get(i) + " ";
			}
			output += "\n";
		}
		System.out.println("done");
		return output;
	}

	public static String coalesceHistograms(int n, int trials) {
		String output = "";
		ArrayList<ArrayList<Integer>> arr = new ArrayList<>();
		for (int i = 3; i <= n; i++) {
			output += i + " ";
			arr.add(coalesceTimeList(i, 1, 0, trials));
			// give status update.
			if ((n > 3 && i % (n / 4) == 0) || i == n) {
				System.out.println(i);
			}
		}
		System.out.println("printing");
		output += "\n";
		for (int i = 0; i < trials; i++) {
			for (int j = 0; j < arr.size(); j++) {
				output += arr.get(j).get(i) + " ";
			}
			output += "\n";
		}
		System.out.println("done");
		return output;
	}

	private static int total(HashMap<Integer, Integer> map) {
		int total = 0;
		for (Entry<Integer, Integer> str : map.entrySet()) {
			total += str.getValue();
		}
		return total;
	}

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
		s.append("\nBinomial Mixture M" + mixtureNToString(n, m, prop));
		s.append("\nConditional Distribution M|N\n" + conditionalsOfAllN(n, m, map));
		s.append("\nConditional Distribution N|M\n" + conditionalsOfAllM(n, m, map));
		return s.toString();
	}

	public static double timeToAbsorbtion(int amount, int n, int m, String state, boolean voter) {
		int[] pars = parse(state);
		int count = 0;
		double avg = 0.0;
		while (count < amount) {
			BipartiteGraph graph = new BipartiteGraph(n, m, voter);
			graph.set(pars[0], pars[1]);
			graph.simulation();
			avg += (double) graph.time() / amount;
			count++;
		}
		return avg;
	}

	public static double proportionOfTimeIn0orN(int n, int m, int steps, int mult, boolean voter) {
		HashMap<String, Integer> map = quasiStationaryDistributionHash(n, m, steps, mult, voter);
		int amount = 0;
		for (Entry<String, Integer> str : map.entrySet()) {
			int[] state = parse(str.getKey());
			if (state[0] == 0 || state[0] == n) {
				amount += str.getValue();
			}
		}
		return ((double) amount / steps) / mult;
	}

	public static String proportionOfTimeIn0orNPlot(int maxn, int maxm, int steps, int mult, boolean voter) {
		HashMap<Integer, HashMap<Integer, Double>> map = new HashMap<>();
		for (int i = 3; i <= maxn; i++) {
			map.put(i, new HashMap<>());
			for (int j = 2; j <= Math.min(i, maxm); j++) {
				map.get(i).put(j, proportionOfTimeIn0orN(i, j, steps, mult, voter));
			}
		}
		String s = "(n,m)   ";
		for (int i = 2; i <= maxm; i++) {
			s += i + "\t";
		}
		for (Entry<Integer, HashMap<Integer, Double>> ent : map.entrySet()) {
			s += "\n" + ent.getKey() + "\t";
			for (Entry<Integer, Double> val : ent.getValue().entrySet()) {
				s += val.getValue() + "\t";
			}
		}
		return s;
	}

	public static void main(String[] args) {
		try {
			PrintWriter out = new PrintWriter("src/voterInvasionModels/Invasion.txt");
			int n = 300;
			int m = 10;
			out.println(allQSDData(n, m, 100000, 10, false));
//			BipartiteGraph g = new BipartiteGraph(200, 1, false);
//			g.simulation();
//			out.println(g.toString());
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
