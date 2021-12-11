package voterInvasionModels;

import java.util.ArrayList;

public class Person {

	private int opinion;
	private int personID;
	private ArrayList<Person> neighbors = new ArrayList<>();

	public Person(int id, int op) {
		opinion = op;
		personID = id;
	}

	public Person(int id) {
		personID = id;
	}

	public void setOpinion(int op) {
		opinion = op;
	}

	public int getOpinion() {
		return opinion;
	}

	public int getID() {
		return personID;
	}

	public boolean addNeighbor(Person p) {
		if (!neighbors.contains(p) && !p.neighbors.contains(this)) {
			neighbors.add(p);
			p.neighbors.add(this);
			return true;
		}
		return false;
	}

	public boolean addNeighbors(ArrayList<Person> people) {
		boolean added = true;
		for (Person p : people) {
			if (!addNeighbor(p)) {
				added = false;
			}
		}
		return added;
	}

	public boolean removeNeighbor(Person p) {
		return neighbors.remove(p);
	}

	public boolean removeNeighbors(ArrayList<Person> people) {
		boolean removed = true;
		for (Person p : people) {
			if (!removeNeighbor(p)) {
				removed = false;
			}
		}
		return removed;
	}

	public ArrayList<Person> neighbors() {
		return neighbors;
	}

	protected void sortneighbors() {
		neighbors.sort((x, y) -> x.getID() - y.getID());
	}

	@Override
	public String toString() {
		return "[" + personID + "] op: " + opinion;
	}

	@Override
	public boolean equals(Object person) {
		if (person == null || !(person instanceof Person)) {
			return false;
		}
		if (personID == ((Person) person).getID()) {
			return true;
		}
		return false;
	}

	@Override
	public int hashCode() {
		return personID;
	}

	public Person clone() {
		return new Person(personID, opinion);
	}
}
