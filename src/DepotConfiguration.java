import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import MDVRP.Node;

public class DepotConfiguration implements Comparable<DepotConfiguration> {
	
	TreeSet<DepotOption> openDepots;
	double costs;
	
	
	public DepotConfiguration copy() {
		DepotConfiguration copy = new DepotConfiguration();
		copy.setCosts(costs);
		for (DepotOption depot : openDepots) {
			copy.add(depot.copy());
		}
		return copy;
	}
	
	
	@Override  
	public int hashCode() {		
		final int prime = 31;
	    int result = 1;
	    for (DepotOption option : openDepots) {
		    result = prime * result + ((option == null) ? 0 : option.hashCode());
	    }
	    return result;
	}
	
	
	public boolean equals(Object obj) {
        return (obj instanceof DepotConfiguration && ((DepotConfiguration) obj).hashCode() == this.hashCode());
    }
	
	
	public DepotConfiguration() {
		openDepots = new TreeSet<DepotOption>();
	}
	
	
	public DepotConfiguration(DepotOption depot) {
		openDepots = new TreeSet<DepotOption>();
		openDepots.add(depot);
	}
	
	
	public void setCosts(double costs) {
		this.costs = costs;
	}
	
	
	public double getCosts() {
		return this.costs;
	}
	
	
	public boolean hasSufficientCapacity(int demand) {
		int sumCapacity = openDepots.stream()
				.mapToInt(DepotOption::getCapacity)
				.sum();	
		return sumCapacity >= demand;
	}
	
	
	public void add(DepotOption depot) {
		openDepots.add(depot);
	}
	
	
	public ArrayList<Node> toMDVRPNodes() {		
		ArrayList<Node> nodes = new ArrayList<Node>();
		
		for (DepotOption depot: openDepots) {
			nodes.add(depot.toMDVRPNode().copy());
		}
		return nodes;
	}
	
	
	public double getOpeningCosts() {
		return openDepots.stream()
				.mapToDouble(DepotOption::getOpeningCost)
				.sum();		
	}

	
	public int compareTo(DepotConfiguration otherConfig) {
		if (Double.compare(costs, otherConfig.costs) == 0)
			return Integer.compare(hashCode(), otherConfig.hashCode());
		else
			return Double.compare(costs, otherConfig.costs);
	}

}
