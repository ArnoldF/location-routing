import java.util.List;
import java.util.concurrent.Callable;

import MDVRP.MDVRPModel;
import MDVRP.Stopwatch;

class SolverThread implements Callable<List<DepotConfiguration>> {
		
		MDVRPModel mdvrp;
		List<DepotConfiguration> configurations;
		int maxIterations;

		public SolverThread(MDVRPModel mdvrp, List<DepotConfiguration> configurations, Integer maxIterations) {
			this.mdvrp = mdvrp;
			this.configurations = configurations;
			this.maxIterations = maxIterations;
		}
		
        @Override
        public List<DepotConfiguration> call() throws Exception {
            try {
            	Stopwatch stopwatch = new Stopwatch();
    			stopwatch.start();
    			//this.setPriority(10);
   			
            	for (DepotConfiguration config : this.configurations) {  		
            		this.mdvrp.setDepots(config.toMDVRPNodes());
            		this.mdvrp.constructStartingSolution();
            		this.mdvrp.checkSolution();

	        		if (maxIterations > 0) {
	        			int maxTime = maxIterations*mdvrp.customers*20;
	        			if (mdvrp.customers >= 1000)
	        				maxTime = maxIterations*mdvrp.customers/2;
	        			this.mdvrp.optimizeRoutes(maxIterations, maxTime, false);
	        		}	        		
	        		double routingCost = mdvrp.computeSolutionCosts();
	        		config.setCosts(routingCost + config.getOpeningCosts());        			 
            	}
            	stopwatch.stop();
    			//System.out.println(stopwatch.elapsedTimeMillis() / 1000.0 + " " +this.configurations.size());
    			return this.configurations;
        		
            }catch (Exception ex){
                ex.printStackTrace();}
            return null;
        }
	}

