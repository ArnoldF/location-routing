import java.io.File;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.StandardOpenOption;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.JsonArray;
import org.json.simple.JsonObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import MDVRP.MDVRPModel;
import MDVRP.Node;

import org.apache.log4j.Logger;


public class LRP {
	
	static final int BASERUNTIME = 60000;	
	static final int MAXCONFIGS = 20000;
	private int maxThreads = 1;
	MDVRPModel mdvrp;
	List<DepotOption> possibleDepots;
	TreeSet<DepotConfiguration> remainingConfigurations;
	int numPossibleDepots;
	int sumDemand;
	int upperBoundOpenDepots;
	private static Logger logger = Logger.getLogger(LRP.class);
	
	
	public LRP(int maxThreads) {
	
		remainingConfigurations = new TreeSet<DepotConfiguration>();	
		this.maxThreads = maxThreads;
	}
	
	
	public void loadInstance(String path) {
		
		try (FileReader reader = new FileReader(path)) {
			JSONParser parser = new JSONParser();
			JSONObject jsonObject = (JSONObject) parser.parse(reader);
			
			mdvrp = new MDVRPModel(jsonObject);		
			if (possibleDepots == null || possibleDepots.isEmpty())
				loadPossibleDepots(jsonObject);
			numPossibleDepots = possibleDepots.size();
			
			sumDemand = Arrays.stream(mdvrp.nodes)
				.mapToInt(n -> n.demand)
				.sum();
			
			logger.info("Instance " + path + " has been loaded." );
				
		} catch (IOException e) {
			logger.error("Instance " + path + " could not be found. Check path.");
	    } catch (ParseException e) {
	    	logger.error("Instance did not match input format.");
	    }	
	}
	
		
	private void loadPossibleDepots(JSONObject jsonObj) throws ParseException {
		
		possibleDepots = new ArrayList<DepotOption>();
		JSONArray depotInfo = (JSONArray) jsonObj.get("depots");     
	    Iterator<Object> iterator = depotInfo.iterator();
	    int index = 0;
	    while (iterator.hasNext()) {
	    	JSONObject nextEntry = (JSONObject) iterator.next();
	    	int capacity = ((Long) nextEntry.get("capacity")).intValue();
	    	double cost = ((Long) nextEntry.get("costs")).intValue();
	    	double x = ((Long) nextEntry.get("x")).intValue();
	    	double y = ((Long) nextEntry.get("y")).intValue();
	    	possibleDepots.add(new DepotOption(index, x, y, cost, capacity));
	    	index++;
	    }
	}
	
	
	public void solve(Integer[] filterStages) {
		
		try {
			upperBoundOpenDepots =  computeUpperBoundOpenDepots();

			System.out.println("Upper bound of open depots has been computed as " + upperBoundOpenDepots + ".");
			
			generateAllDepotConfigurations(upperBoundOpenDepots);
			if (remainingConfigurations.size() > MAXCONFIGS) {
				logger.info("Apply heuristic construction.");
				heuristicConstruction();
			}	
			System.out.println("All " + remainingConfigurations.size() +" relevant depot configurations have been build.");
			
			applyFilters(filterStages);
			logger.info("Optimation completed.");
		}
		catch (Exception e) {
			logger.error("An error occured during optimisation.");
		}
	}
	
	
	private int computeUpperBoundOpenDepots() {
		
		List<Integer> sortedCapacity = possibleDepots.stream()
				.map(d -> d.getCapacity())
				.sorted()
				.collect(Collectors.toList());
				
		//determine minimal number of open depots, given limited depot capacity
		int sumCapacity = 0;
		int minOpenDepots = 0;
		while (minOpenDepots < numPossibleDepots && sumCapacity < sumDemand){
			sumCapacity += sortedCapacity.get(minOpenDepots);
			minOpenDepots++;
		}
		
		//determine most central depot
		double shortestDistance = Double.MAX_VALUE;
		DepotOption centralDepot = null;
		for (DepotOption d: possibleDepots) {
			double sumDistance = 0;
			for (Node n: mdvrp.nodes) {
				sumDistance += (float) Math.sqrt( Math.pow( n.x - d.getX() , 2) + Math.pow( n.y - d.getY() , 2) );				
			}
			if (sumDistance < shortestDistance) {
				centralDepot = d;
				shortestDistance = sumDistance;
			}
		}
		
		//estimate routing costs from central depot
		mdvrp.setDepots(new DepotConfiguration(centralDepot).toMDVRPNodes());
		mdvrp.depotInventory[0][0] = Integer.MAX_VALUE;
		double r1 = mdvrp.constructStartingSolution();
		
		List<Double> sortedOpeningCosts = possibleDepots.stream()
				.map(d -> d.getOpeningCost())
				.sorted()
				.collect(Collectors.toList());	
		int minNumRoutes = (int)Math.ceil((double)sumDemand / mdvrp.capacityLimitVehicle);
		
		int upperBound = 2;
		double costReduction = r1 * (estimatedRoutingReduction(upperBound, minNumRoutes) - estimatedRoutingReduction(upperBound-1, minNumRoutes));
		while(upperBound < numPossibleDepots && costReduction > sortedOpeningCosts.get(upperBound-1)) {
			upperBound++;
			costReduction = r1 * (estimatedRoutingReduction(upperBound, minNumRoutes) - estimatedRoutingReduction(upperBound-1, minNumRoutes));
		}
		mdvrp.setDepots(new ArrayList<Node>());
		return Math.min(Math.max(upperBound, minOpenDepots), 200);
	}
	
	
	private double estimatedRoutingReduction(int M, int minNumRoutes) {
		
		double val = Math.pow(Math.log((double)M), 0.58) * 0.27 * Math.pow(minNumRoutes / 10.0, 1.0 / 3);
		return val;
	}
	
	
	private void generateAllDepotConfigurations(int upperBoundOpenDepots) {
		
		remainingConfigurations = new TreeSet<DepotConfiguration>();
		
		generationLoop:
		for (int numOpenDepots = 1; numOpenDepots <= upperBoundOpenDepots; numOpenDepots++) {
			int[] current_depot_option = new int[numOpenDepots];
			for (int i=0; i<numOpenDepots; i++)
				current_depot_option[i] = i;
			while (true) {			
				DepotConfiguration currentDepotConfig = new DepotConfiguration();
				for (int i=0; i<numOpenDepots; i++) {
					currentDepotConfig.add(possibleDepots.get(current_depot_option[i]));	
				}
				if (currentDepotConfig.hasSufficientCapacity(sumDemand)) {
					remainingConfigurations.add(currentDepotConfig);
				}
				if (remainingConfigurations.size() > MAXCONFIGS)
					break generationLoop;
	
				int posToChange = numOpenDepots - 1;
				while ( posToChange >= 0 && current_depot_option[posToChange] == numPossibleDepots - (numOpenDepots - posToChange))
					posToChange--;
				if (posToChange == -1)
					break;
				current_depot_option[ posToChange ]++;
				int resetValue = current_depot_option[ posToChange ];
				for (int i=posToChange+1; i<numOpenDepots; i++)
					current_depot_option[ i ] = resetValue + i - posToChange;				
			}
		}
	}		
	
	
	private void applyFilters(Integer[] filterStages) {
			
		try {
			evaluateDepotConfigurations(0);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}//TODO pass threshold as parameter to abort CW early
		System.out.println("Best objective value: " + remainingConfigurations.first().costs);
		
		for (int numCap: filterStages) {
			filterRemainingOptions(numCap);
			double timeLimit = BASERUNTIME / (100 * Math.max(1, filterStages.length-1) * numCap);
			int maxIterations = (int)Math.ceil(timeLimit * 10);
			if (numCap <= maxThreads) {
				maxIterations = 3000;
				if (mdvrp.customers >= 1000)
					maxIterations = 100;
			}

			try {
				evaluateDepotConfigurations(maxIterations);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			System.out.println("The best " + numCap + " configurations have been evaluated. Best objective value: " + remainingConfigurations.first().costs);
			if (numCap <= maxThreads)
				break;
		}		
	}
	

	private void filterRemainingOptions(int maxOptions)	{
		
		while (remainingConfigurations.size() > maxOptions) {
			remainingConfigurations.remove(remainingConfigurations.last());
		}
	}


	private void evaluateDepotConfigurations(int maxIterations) throws InterruptedException {
		
		//divide remainingDepotConfigurations into equal parts
		int configurationsPerThread = (int) Math.ceil((double)remainingConfigurations.size() / maxThreads);
		ExecutorService executor = Executors.newFixedThreadPool(maxThreads);
		List<Future<List<DepotConfiguration>>> resultList = new ArrayList<>();
		
		ArrayList<ArrayList<DepotConfiguration>> configurations = new ArrayList<ArrayList<DepotConfiguration>>();
		configurations.add(new ArrayList<DepotConfiguration>());
		int countAssignedConfigs = 0;
		int threadIndex = 0;
		for (DepotConfiguration config : remainingConfigurations) {
			configurations.get(threadIndex).add(config.copy());
			countAssignedConfigs++;
			if (countAssignedConfigs >= configurationsPerThread) {
				countAssignedConfigs = 0;
				threadIndex ++;
				configurations.add(new ArrayList<DepotConfiguration>());
			}
		}
		for (int threadIdx= 0; threadIdx < maxThreads; threadIdx++) {
			SolverThread newthread = new SolverThread(mdvrp.copy(), configurations.get(threadIdx), maxIterations);		
			Future<List<DepotConfiguration>> result1 = executor.submit(newthread);
			resultList.add(result1);
		}
                     
		remainingConfigurations.clear();
		
		executor.shutdown();
	    try {
	        if (!executor.awaitTermination(600, TimeUnit.SECONDS)) {
	        	executor.shutdownNow();
	        }
	    } catch (InterruptedException ex) {
	    	executor.shutdownNow();
	        Thread.currentThread().interrupt();
	    }
				
		for (Future<List<DepotConfiguration>> result: resultList) {
			try {
				remainingConfigurations.addAll(result.get());
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
	}
		
	
	public void writeToFile() {
		
		JSONObject output = new JSONObject();
		output.put("Objective", remainingConfigurations.first().costs);
		output.put("NumOpenDepots", remainingConfigurations.first().openDepots.size());

	    JSONArray openDepots = new JSONArray();
	    for (DepotOption depot : remainingConfigurations.first().openDepots) {
	    	openDepots.add(depot.getIndex());
	    }
	    output.put("OpenDepots", openDepots);
	    
	    Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat("HH_mm_ss");
        Random rand = new Random();
        int randomNum = rand.nextInt(9999);
        String path = "." + File.separator +"output" + File.separator + "Result_"+sdf.format(cal.getTime()) + String.valueOf(randomNum) + ".sol";

	    try (FileWriter file = new FileWriter(path, true)) {
	        file.write(output.toJSONString());
	        logger.info("Written solution.");
	        logger.info(path);
	    } catch (IOException e) {
	    	logger.error("Output File could not be written. Check whether the output folder exists.");
	    }
	}
	
	
	private void heuristicConstruction() {

		remainingConfigurations.clear();
		for (int r=1; r<=Math.max(16, this.upperBoundOpenDepots); r++) {
			HashSet<DepotConfiguration> newConfigs = HeuristicConstruction.createConfigs(this, r, MAXCONFIGS);
			for (DepotConfiguration config : newConfigs)
				remainingConfigurations.add(config);
			if (remainingConfigurations.size() >= this.MAXCONFIGS)
				break;
		}
	}
	

}
