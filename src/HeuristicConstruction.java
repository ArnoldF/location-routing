import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import MDVRP.Node;

public class HeuristicConstruction {
	
	static private double minY;
	static private double minX;
	static private double maxY;
	static private double maxX;
	static private Map<DepotOption, Double> depotGrossCosts;
	static private ArrayList<Region> regions;
	static private int upperBoundOpenDepots;
	static private int sumDemand;
	
	
	public static class Region {
		Double x;
		Double y;
		Double width;
		Double height;
		
		public Region(Double x, Double y, Double width, Double height) {
			this.x = x;
			this.y = y;
			this.width = width;
			this.height = height;
		}
		
		public boolean contains(Double pointX, Double pointY) {
			if (pointX >= x - (width / 2.0) && pointX < x + (width / 2.0))
				if (pointY >= y - (height / 2.0) && pointY < y + (height / 2.0))
					return true;
			return false;
		}
	}
	
	
	public static HashSet<DepotConfiguration> createConfigs(LRP lrp, int numRegions, int maxConfigurations) {
		
		upperBoundOpenDepots = lrp.upperBoundOpenDepots;
		sumDemand = lrp.sumDemand;
		computeDepotGrossCosts(lrp);
		findEdgeCoordinates(lrp);		
		computeRegions(numRegions);
		LinkedList<Set<Set<DepotOption>>> depotSubConfigurations = createSubConfigurations(lrp);	
		return createConfigurations(depotSubConfigurations);		
	}
	
	
	private static void findEdgeCoordinates(LRP lrp) {
		
		minY =	Double.MAX_VALUE;
		minX =	Double.MAX_VALUE;
		maxY = -Double.MAX_VALUE;
		maxX = -Double.MAX_VALUE;
		for (DepotOption config : lrp.possibleDepots) {
			minX = Math.min(minX, config.getX());
			minY = Math.min(minY, config.getY());
			maxX = Math.max(maxX, config.getX());
			maxY = Math.max(maxY, config.getY());
		}
		for(Node n: lrp.mdvrp.nodes){
			minX = Math.min(minX, n.x);
			minY = Math.min(minY, n.y);
			maxX = Math.max(maxX, n.x);
			maxY = Math.max(maxY, n.y);
		}
		maxY = maxY*1.001;
		maxX = maxX*1.001;
	}
	
	
	private static void computeDepotGrossCosts(LRP lrp) {
		
		depotGrossCosts = new TreeMap<DepotOption, Double>();
		double numRoutesEstimate = Math.ceil( lrp.sumDemand / lrp.mdvrp.capacityLimitVehicle );
		int avgToursPerDepot = (int) Math.ceil( numRoutesEstimate / lrp.upperBoundOpenDepots );
	
		for (DepotOption depot : lrp.possibleDepots) {

			ArrayList<Double> distancesToDepotSorted = new ArrayList<Double>();
			for(Node customer : lrp.mdvrp.nodes) {
				//TODO
				//distancesToDepotSorted.add( (double) Math.ceil(100.0 * Math.sqrt(Math.pow(depot.getX() - customer.x , 2) + Math.pow(depot.getY() - customer.y , 2) )));	
				distancesToDepotSorted.add( (double)Math.sqrt(Math.pow(depot.getX() - customer.x , 2) + Math.pow(depot.getY() - customer.y , 2) ));
			}	
			Collections.sort(distancesToDepotSorted);
						
			double grossCosts = depot.getOpeningCost();
			for (int distanceIndex=0; distanceIndex<Math.min(avgToursPerDepot*2, lrp.mdvrp.customers); distanceIndex++)
				grossCosts += distancesToDepotSorted.get(distanceIndex);
			grossCosts = grossCosts / (float)depot.getCapacity();
			depotGrossCosts.put(depot, grossCosts);
		}
		depotGrossCosts = (TreeMap<DepotOption, Double>) sortByValues(depotGrossCosts);
	}
	
	
	private static void computeRegions(int numRegions) {
	
		regions = new ArrayList<Region>();
		int cols = (int) Math.ceil( Math.sqrt(numRegions) );
		int rows = (int) Math.ceil((double)numRegions / (double)cols);
		int numShortRows = (cols * rows) - numRegions;

		double widthTotal = maxX - minX;
		double heightTotal = maxY - minY;
		double widthRegion 	= widthTotal / cols;
		double heightRegion	= heightTotal / rows;
					
		double y	= minY + heightRegion/2;
		boolean prevRowIsShort = false;
		if ((double)numShortRows < (double)rows/2) {
			prevRowIsShort = true;
		}		
					
		for (int rowCounter=0; rowCounter<rows; rowCounter++) {
			int colsInRow = cols;
			double x = minX + widthRegion/2;
			if (numShortRows > 0 && !prevRowIsShort) {
				//create short row with one less region
				x = minX + widthRegion;
				colsInRow = cols-1;
				numShortRows--;
				prevRowIsShort = true;
			}
			else {				
				prevRowIsShort = false;
			}											
			for (int countCols=0; countCols<colsInRow; countCols++) {
				regions.add( new Region(x, y, widthRegion, heightRegion));
				x += widthRegion;
			}
			y += heightRegion;
		}
	}
	
	
	private static LinkedList<Set<Set<DepotOption>>> createSubConfigurations(LRP lrp) {
				
		ArrayList<DepotOption> depotsSortedByCosts = new ArrayList<DepotOption>();
		for (DepotOption depot : depotGrossCosts.keySet()) {
			depotsSortedByCosts.add(depot);
		}			
		LinkedList<Set<Set<DepotOption>>> depotSubConfigurations = new LinkedList<Set<Set<DepotOption>>>();
		
		//find depots per region, sorted by grossCosts ascending
		for (Region region : regions) {
			ArrayList<DepotOption> depotsInRegion = new ArrayList<DepotOption>();
			for (DepotOption depot : depotsSortedByCosts) {
				if (region.contains(depot.getX(), depot.getY())) {
					depotsInRegion.add(depot);
					if (lrp.mdvrp.customers >= 1000)
						break;
				}
			}
			Set<Set<DepotOption>> depotSubConfigurationsPerRegion = new LinkedHashSet<Set<DepotOption>>();
			Set<DepotOption> subConfiguration = new LinkedHashSet <DepotOption>();
			if (lrp.mdvrp.customers < 1000) //add empty set
				depotSubConfigurationsPerRegion.add(new LinkedHashSet <DepotOption>(subConfiguration));
			for (DepotOption depot : depotsInRegion) {
				subConfiguration.add(depot);
				depotSubConfigurationsPerRegion.add(new LinkedHashSet <DepotOption>(subConfiguration));				
			}
			depotSubConfigurations.add(depotSubConfigurationsPerRegion);
		}
		return depotSubConfigurations;
	}
	
	
	private static HashSet<DepotConfiguration> createConfigurations(LinkedList<Set<Set<DepotOption>>> depotSubConfigurations) {

		ListIterator<Set<Set<DepotOption>>> iterator = depotSubConfigurations.listIterator(); 
		HashSet<DepotConfiguration> configurations = new HashSet<DepotConfiguration>();
		combineSubConfigurations( depotSubConfigurations, iterator, new LinkedHashSet<Set<DepotOption>>(), configurations );
		return configurations;	
	}
	
	
    private static void combineSubConfigurations( LinkedList<Set<Set<DepotOption> >> subConfigurations, ListIterator<Set<Set<DepotOption>>> iter, Set<Set<DepotOption>> cur, Set<DepotConfiguration> configurations ) {

        if (!iter.hasNext()) {
        	DepotConfiguration config = new DepotConfiguration();

            for (Set<DepotOption> set : cur) {
            	for (DepotOption depot : set) {
            		config.add(depot);
            	}
            }
            if (config.openDepots.size() <= upperBoundOpenDepots && config.hasSufficientCapacity(sumDemand) ) {
            	configurations.add( config );
            }
        } else {
        	Iterator<Set<Set<DepotOption>>> previous = iter;
        	Set<Set<DepotOption>> set = iter.next();
        	
        	if (set.isEmpty())
        		combineSubConfigurations( subConfigurations, iter, cur, configurations );
        	else{
	            for( Set<DepotOption> value : set ) {
	                cur.add( value );
	                combineSubConfigurations( subConfigurations, iter, cur, configurations );
	                cur.remove( value );
	            }
        	}
            iter.previous();
        }
    }
    
    
	private static <K, V extends Comparable<V>> Map<K, V> sortByValues(final Map<K, V> map) {
	    Comparator<K> valueComparator = new Comparator<K>() {
	    	public int compare(K k1, K k2) {
	    		int compare = map.get(k1).compareTo(map.get(k2));
	    		if (compare == 0) 
	    			return 1;
	    		else 
	    			return compare;
	    	}
	    };
	 
	    Map<K, V> sortedByValues = new TreeMap<K, V>(valueComparator);
	    sortedByValues.putAll(map);
	    return sortedByValues;
  }
}
