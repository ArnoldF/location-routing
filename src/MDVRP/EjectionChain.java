package MDVRP;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class EjectionChain {

	ArrayList<ArrayList<Relocation>> incomplete_moves;
	ArrayList<ArrayList<Relocation>> candidate_moves;
	ArrayList<Relocation> relocations_all;
	MDVRPModel model;
	
	
	ArrayList<Integer> relocations_to_be_extended = new ArrayList<Integer>();
	
	private class Relocation
	{
		
		public double saving_p;
		public double saving;
		public Node node;
		public int fromRoute;
		public int toRoute;
		public int newpos;
		public int chainPos;
		public int buffer;
		public double bufferTour;
		public double oldDetour;
		public double newDetour;
		public Relocation predecessor;
		public ArrayList<Integer> forbiddenNodes;
		public ArrayList<Integer> relocateDestinations;
			
		public Relocation(double v, double vT, Node n, int f, int t, int p, int c, int b, double bT, double oD, double nD, Relocation pre, ArrayList<Integer> forbidden, ArrayList<Integer> destinations)
		{
			saving_p 		= v;
			saving		= vT;
			node 		= n;
			fromRoute 	= f;
			toRoute		= t;
			newpos		= p;
			chainPos	= c;
			buffer		= b;
			bufferTour  = bT;
			oldDetour 	= oD;
			newDetour	= nD;
			predecessor = pre;
			forbiddenNodes = (ArrayList<Integer>) forbidden.clone();
			relocateDestinations = (ArrayList<Integer>) destinations.clone();
		}
		
		public Relocation() {
			// TODO Auto-generated constructor stub
			forbiddenNodes = new ArrayList<Integer>();
			relocateDestinations = new ArrayList<Integer>();
		}

		public void clear()
		{
			node = null;
			forbiddenNodes.clear();
		}
		
		public Relocation copy()
		{
			Relocation copy = new Relocation(this.saving_p, this.saving, this.node, this.fromRoute, this.toRoute, this.newpos, this.chainPos, this.buffer, this.bufferTour, this.oldDetour, this.newDetour, this.predecessor, this.forbiddenNodes, this.relocateDestinations );
			return copy;
		}
	}//class
	
	private static class MapUtil {
	    public static <K, V extends Comparable<? super V>> Map<K, V> 
	        sortByValue(Map<K, V> map) {
	        List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>(map.entrySet());
	        Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
	            public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2) {
	                return (o1.getValue()).compareTo( o2.getValue() );
	            }
	        });

	        Map<K, V> result = new LinkedHashMap<K, V>();
	        for (Map.Entry<K, V> entry : list) {
	            result.put(entry.getKey(), entry.getValue());
	        }
	        return result;
	    }
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//FUNCTION : localSearchInterRoute
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	public void ejectionChain( MDVRPModel modelOriginal, int max_relocations, ArrayList<Integer> nodesToStartChains )
	{		
		//INITIALIZATION
		model = modelOriginal;
		double maxWorsening = -10 * model.baselineUnit;
		incomplete_moves = new ArrayList<ArrayList<Relocation>>();
		candidate_moves = new ArrayList<ArrayList<Relocation>>();
		relocations_all = new ArrayList<Relocation>();

		ArrayList<Integer> toursToContinueChainsFrom = new ArrayList<Integer>();
		
		//----------------------------------------
		// FIND STARTING RELOCATIONS FOR RELOCATION CHAINS
		//----------------------------------------
		for (int p=0; p<nodesToStartChains.size(); p++)
		{			
			 Node starting_node = model.nodes[ nodesToStartChains.get(p) ];
			 //sort the nearest nodes of starting_node according to insertion costs
			 starting_node.insert_next_to = (HashMap<Integer, InsertionCost>) MapUtil.sortByValue(starting_node.insert_next_to);
			Iterator insertions = starting_node.insert_next_to.entrySet().iterator();		
			ArrayList<Integer> insertion_into = new ArrayList<Integer>();
			
			//allow to relocate the starting_node into an empty route (only for MDVRPs)
			if (model.depots > 1)
			if ( starting_node.detour_p > 2 * model.distance_to_nearest_depot[starting_node.index].costs + model.cost_per_Route  && starting_node.detour - 2 * model.distance_to_nearest_depot[starting_node.index].costs > maxWorsening )
			{
				int targetDepot = model.distance_to_nearest_depot[starting_node.index].node1;
				//find empty route of this depot
				for (int emptyRoute=0; emptyRoute < model.routes.size(); emptyRoute++)
				if ( model.routes.get(emptyRoute).length == 0 && model.routes.get(emptyRoute).depot == targetDepot )
				{
					extendRelocationChain( starting_node, starting_node.routeIndex, emptyRoute, 1, starting_node.detour_p - 2 * model.distance_to_nearest_depot[starting_node.index].costs - model.cost_per_Route, starting_node.detour - 2 * model.distance_to_nearest_depot[starting_node.index].costs - model.cost_per_Route, starting_node.detour, 2 * model.distance_to_nearest_depot[starting_node.index].costs, null );		
				}
			}
			
			//try to insert starting_node next to its nearest nodes
		    while (insertions.hasNext())
		    {
				relocations_all.clear();
				relocations_to_be_extended.clear();
		        Map.Entry entry = (Map.Entry) insertions.next();
				InsertionCost insertion = (InsertionCost)entry.getValue();
				
				//compute the costs of a relocation next to the near node
				double saving_p = starting_node.detour_p - insertion.costs_p;
				double saving 	= starting_node.detour - insertion.costs;
				//LRP
                if (model.routes.get(starting_node.routeIndex).length == 1)
                {
                	saving_p += model.cost_per_Route;
                	saving 	 += model.cost_per_Route;
                }


				model.ev3++;
				if ( saving_p >= 0 && saving >= maxWorsening && !insertion_into.contains(insertion.route))
				{
					//add a relocation as potential start to a chain
				insertion_into.add(insertion.route);
				int newRoute  			= insertion.route;
				int newPos  			= insertion.position;
				extendRelocationChain( starting_node, starting_node.routeIndex, newRoute, newPos, saving_p, saving, starting_node.detour, insertion.costs, null );
				
				//-------------------------------------------------------------------------
				//TRY TO EXTEND RELOCATION CHAINS
				//-------------------------------------------------------------------------
				loopRoutes:
				for (int level = 1; level < max_relocations; level++)
				{
					//relocations_to_be_extended is updated in function extendEjectionChain()
					ArrayList<Integer> relocations_to_be_extended_copy = (ArrayList<Integer>) relocations_to_be_extended.clone();
					relocations_to_be_extended.clear();
					
					//try to extend all incomplete relocation chains
					for (int nextIndex : relocations_to_be_extended_copy)						
					{
						Relocation nextOption = relocations_all.get( nextIndex );
						Route fromRoute = model.routes.get(nextOption.toRoute);
						//try to relocate each node in the destination route
						for(Node nextNode : fromRoute.nodes)
						if (nextNode.index < model.customers)
						{
							//if we would remove this node from the route, is the route feasible?
							int bufferCapacity =  nextOption.buffer + nextNode.demand;
							double bufferTourLength = nextOption.bufferTour + nextNode.detour;
							if ( bufferCapacity >= 0 && bufferTourLength >= 0 )
							{
								//try to relocate this node next to its nearest nodes
								Iterator insertions_continued = nextNode.insert_next_to.entrySet().iterator();		
								ArrayList<Integer> insertion_into_continued = new ArrayList<Integer>();
							    while (insertions_continued.hasNext())
							    {
							    	entry = (Map.Entry)insertions_continued.next();
									insertion = (InsertionCost)entry.getValue();						    	
									model.ev3++;
									
									//check whether the sum of all relocations in the chain results in a saving
									saving_p		= nextOption.saving_p  + nextNode.detour_p - insertion.costs_p;
									saving 			= nextOption.saving + nextNode.detour - insertion.costs;
									newRoute  		= insertion.route;
									newPos  		= insertion.position;
									//if a node next to the insertion position of the previous relocation is relocated, costs have to be corrected
									if (nextOption.newpos == nextNode.routePosition)
									{
										saving_p += model.getDistance2( nextNode.index, nextOption.node.index ) + model.getDistance2( nextNode.neighbours[0].index, nextNode.neighbours[1].index ) -  model.getDistance2( nextNode.index, nextNode.neighbours[0].index ) - model.getDistance2( nextOption.node.index, nextNode.neighbours[1].index );
										saving += model.getDistance( nextNode.index, nextOption.node.index ) + model.getDistance( nextNode.neighbours[0].index, nextNode.neighbours[1].index ) - model.getDistance( nextNode.index, nextNode.neighbours[0].index ) - model.getDistance( nextOption.node.index, nextNode.neighbours[1].index );
									}
									else if (nextOption.newpos == nextNode.routePosition + 1)
									{
										saving_p +=  model.getDistance2( nextNode.index, nextOption.node.index ) + model.getDistance2( nextNode.neighbours[0].index, nextNode.neighbours[1].index ) - model.getDistance2( nextNode.index, nextNode.neighbours[1].index ) - model.getDistance2( nextOption.node.index, nextNode.neighbours[0].index );
										saving += model.getDistance( nextNode.index, nextOption.node.index ) + + model.getDistance( nextNode.neighbours[0].index, nextNode.neighbours[1].index ) -  model.getDistance( nextNode.index, nextNode.neighbours[1].index ) - model.getDistance( nextOption.node.index, nextNode.neighbours[0].index );
									}
									
									if ( saving_p >= 0 && saving >= maxWorsening && !insertion_into_continued.contains(newRoute))
									{										
										if (  !nextOption.forbiddenNodes.contains(nextNode.index) &&
												  !nextOption.relocateDestinations.contains(1000*newRoute+newPos) )
										{
											insertion_into_continued.add(newRoute);
											extendRelocationChain( nextNode, nextNode.routeIndex, newRoute, newPos, saving_p, saving, nextNode.detour, insertion.costs, nextOption );
										}
									}//if
								}//while
								
							}//constraints
						}
					}//each node in destination route			    	
				}//length of chain			
			}
			}//insertion positions
			

		}//initial node
		
		//-------------------------------------------------------------------------
		//EXECUTE RELOCATION CHAIN WITH MOST SAVING
		//-------------------------------------------------------------------------
		int chainToBeExecuted = getBestChain();
		while (chainToBeExecuted > -1 )
		{
			//execute the chain
			Relocation executedChain = candidate_moves.get(chainToBeExecuted).get(0);
			executeRelocationChain( chainToBeExecuted );

			//remove all candidate moves that contain forbidden nodes or that now violate constraints
			for (int i=0; i<candidate_moves.size(); i++)
			{
				ArrayList<Relocation> chain = candidate_moves.get(i);
				boolean negativeBuffer = false;
				//inventory constraints
				if (!model.inventoryOff)
				{
					int[][] changedDemand = new int[model.depots][model.products];
					for (int j=0; j<chain.size(); j++)
						for (int p=0; p < model.products; p++)
						{
							changedDemand[model.routes.get(chain.get(j).fromRoute).depot][p] -= chain.get(j).node.demandPerProduct[p];
							changedDemand[model.routes.get(chain.get(j).toRoute).depot][p]   += chain.get(j).node.demandPerProduct[p];
						}
					for (int d=0; d < model.depots; d++)
						if ( !model.hasSufficientInventory(d, changedDemand[d]) )
							negativeBuffer = true;
				}
				//capacity and tour length constraints
				if (chain.get(0).buffer < 0 || chain.get(0).bufferTour < 0)
					negativeBuffer = true;
				else
					for (int j=1; j<chain.size(); j++)
						if (chain.get(j).buffer + chain.get(j-1).node.demand < 0 || chain.get(j).bufferTour + chain.get(j-1).oldDetour < 0 )
						{	negativeBuffer = true; break;}
				if ( negativeBuffer || !Collections.disjoint(chain.get(0).forbiddenNodes, executedChain.forbiddenNodes))
					candidate_moves.get(i).clear();
			}
			candidate_moves.get(chainToBeExecuted).clear();
			Iterator<ArrayList<Relocation>> it = candidate_moves.iterator();
			while (it.hasNext()) 
				if (it.next().isEmpty())
					it.remove();
			
			//find the next best relocation chain
			 chainToBeExecuted = getBestChain();
		}
		incomplete_moves.clear();
		candidate_moves.clear();
	}
	
	
	private Integer getBestChain()
	{
		int indexOfBestChain = -1;
		double minSave = -999999;
		for (int i=0; i<candidate_moves.size(); i++)
		{
			double value = candidate_moves.get(i).get(0).saving_p;
			if ( value > minSave)
			{					
				minSave = value;
				indexOfBestChain = i;
			}
		}
		return indexOfBestChain;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//FUNCTION : execute Injection chain
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	private void executeRelocationChain(int chainIndex)
	{
		model.ex3++;
		//System.out.println( "OUT " + candidate_moves.get(chainIndex).size() +"-RC");
		
		Relocation currentRelocation;
		double prev1 = model.computeSolutionCosts();
		double savePredicted = candidate_moves.get(chainIndex).get(0).saving;
		//execute each relocation of the chain, starting with the last one
		/*if (candidate_moves.get(chainIndex).get(0).toRoute == candidate_moves.get(chainIndex).get( candidate_moves.get(chainIndex).size()-1 ).fromRoute )
			System.out.println("   Y " + candidate_moves.get(chainIndex).size());
		else
			System.out.println("N  ");*/
		for (int chainPos=0; chainPos < candidate_moves.get(chainIndex).size(); chainPos++)
		{
			currentRelocation = candidate_moves.get(chainIndex).get(chainPos);
			int from = currentRelocation.fromRoute;
			Route fromRoute = model.routes.get(from);
			int to   = currentRelocation.toRoute;
			Route toRoute = model.routes.get(to);
			
			//System.out.println( from + " " + currentRelocation.node.routePosition +" "+ to + " " + currentRelocation.newpos + " " + currentRelocation.node.index);
			
			//remove the node from its previous route
			fromRoute.nodes.remove(currentRelocation.node.routePosition);
			fromRoute.length--;
			fromRoute.cost -= currentRelocation.oldDetour;
			fromRoute.load -= currentRelocation.node.demand;
				
			//insert the node in the new route
			toRoute.nodes.add(currentRelocation.newpos, currentRelocation.node);
			toRoute.length++;
			toRoute.cost += currentRelocation.newDetour;
			toRoute.load += currentRelocation.node.demand;
			currentRelocation.node.routeIndex = to;
			
			//update the available inventory, if routes are from different depots
			currentRelocation.node.depot = model.routes.get(to).depot;
			if (!model.inventoryOff && fromRoute.depot != toRoute.depot)
				for (int i=0; i < model.products; i++)
				{
					model.remainingInventory[fromRoute.depot][i] += currentRelocation.node.demandPerProduct[i];
					model.remainingInventory[toRoute.depot][i] -= currentRelocation.node.demandPerProduct[i];
				}
			
			//if we added a new route, increase the penalties
			if (chainPos == 0 && model.routes.get(to).length == 1)
				model.setPenalties(currentRelocation.node.index, model.routes.get(to).nodes.get(0).index, model.getPenalties(currentRelocation.node.index , model.routes.get(from).nodes.get(0).index));
			
			//update neighbours of involved nodes
			currentRelocation.node.neighbours[0].neighbours[1] = currentRelocation.node.neighbours[1];
			currentRelocation.node.neighbours[1].neighbours[0] = currentRelocation.node.neighbours[0];
			currentRelocation.node.neighbours[0] = model.routes.get(to).nodes.get(currentRelocation.newpos-1);
			currentRelocation.node.neighbours[1] = model.routes.get(to).nodes.get(currentRelocation.newpos+1);
			model.routes.get(to).nodes.get(currentRelocation.newpos-1).neighbours[1] = model.nodes[currentRelocation.node.index];
			model.routes.get(to).nodes.get(currentRelocation.newpos+1).neighbours[0] = model.nodes[currentRelocation.node.index];
	
			//post-processing: update the involved routes later
			if (!model.routes_to_be_optimised.contains( to ) )
				model.routes_to_be_optimised.add( to );
			if (!model.routes_to_be_optimised.contains( from ) )
				model.routes_to_be_optimised.add( from );
			if (!model.routes_changed_during_perturbation.contains( to ) )
				model.routes_changed_during_perturbation.add( to );
			if (!model.routes_changed_during_perturbation.contains( from ) )
				model.routes_changed_during_perturbation.add( from );
			
			if (!model.routes_update_RC_variables.contains(to)) 
				model.routes_update_RC_variables.add( to );
			if (!model.routes_update_RC_variables.contains(from)) 
				model.routes_update_RC_variables.add( from );
			
			if (!model.routes_update_edge_values.contains(to)) 
				model.routes_update_edge_values.add( to );
			if (!model.routes_update_edge_values.contains(from)) 
				model.routes_update_edge_values.add( from );
			
			if (!model.penalties_off && !model.heuristicStart)
				model.moves_RC++;

			relocationPostProcessing(from, to, currentRelocation.node.routePosition, currentRelocation.newpos, currentRelocation.node.demand, currentRelocation.oldDetour, currentRelocation.newDetour);
		}	

		//model.checkSolution();
		
		double after = model.computeSolutionCosts();
		//double after2 = this.computeSolutionCosts2();
		//if (prev2 < after2-0.01)
		if ( prev1 - after > savePredicted + 0.1 || prev1 - after < savePredicted - 0.1)//0.0001
			;//System.out.println("ERROR - savings in EC incorrect " + countChains);
		//System.out.println((after - prev1) + " "+ savePredicted);///model.baselineUnit);				
	}
	
	private void relocationPostProcessing(int from, int to, int oldRoutePos, int newRoutePos, int demand, double oldDetour, double newDetour)
	{	
		//update the route position of all customers
		for(int i=1; i < model.routes.get(to).length+1; i++)
			model.routes.get(to).nodes.get(i).routePosition = i;
		for(int i=1; i < model.routes.get(from).length+1; i++)
			model.routes.get(from).nodes.get(i).routePosition = i;
		
		//update the insertion positions of all candidate moves
		for (int i=0; i < candidate_moves.size(); i++)
			for (int j=0; j < candidate_moves.get(i).size(); j++)
				if ( candidate_moves.get(i).get(j).toRoute == from && candidate_moves.get(i).get(j).newpos > oldRoutePos)
					candidate_moves.get(i).get(j).newpos--;
				else if ( candidate_moves.get(i).get(j).toRoute == to && candidate_moves.get(i).get(j).newpos > newRoutePos - 1)
					candidate_moves.get(i).get(j).newpos++;
		
		// adapt the capacity buffer
		for (int i=0; i < candidate_moves.size(); i++)
			for (int j=0; j < candidate_moves.get(i).size(); j++)
				if ( candidate_moves.get(i).get(j).toRoute == from)
					{candidate_moves.get(i).get(j).buffer += demand; candidate_moves.get(i).get(j).bufferTour += oldDetour;}
				else if ( candidate_moves.get(i).get(j).toRoute == to)
					{candidate_moves.get(i).get(j).buffer -= demand; candidate_moves.get(i).get(j).bufferTour -= newDetour;}
									
	}

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//FUNCTION : add a relocation to a relocation chain
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	private void extendRelocationChain ( Node relocate, int oldRoute, int newRoute, int insertionPos, double relocationSaving_p, double relocationSaving, double oldDetour, double newDetour, Relocation bestRelocationOption)
	{	
			//check capacity and route length constraints of the destination route
			int bufferCapacity = model.capacityLimitVehicle - (model.routes.get(newRoute).load + relocate.demand);
			double bufferTour = model.maxTourLength - (model.routes.get(newRoute).cost + newDetour);

			//a preceding relocation might affect these constraints checks
			Relocation temp = bestRelocationOption;
			while ( temp != null)
			{
				if( temp.fromRoute == newRoute)
					{bufferCapacity += temp.node.demand; bufferTour+=temp.oldDetour;}
				else if ( temp.toRoute == newRoute )
					{bufferCapacity -= temp.node.demand; bufferTour-=temp.newDetour;}
				temp = temp.predecessor;
			}
			
			//check inventory constraints
			boolean sufficientInv = true;
			if (!model.inventoryOff)
			{
				int[][] InvChanges = new int[model.depots][model.products];
				for (int p=0; p < model.products; p++)
				{
					InvChanges[model.routes.get(oldRoute).depot][p] -= relocate.demandPerProduct[p];
					InvChanges[model.routes.get(newRoute).depot][p] += relocate.demandPerProduct[p];
				}
				temp = bestRelocationOption;
				while ( temp != null)
				{
					for (int p=0; p < model.products; p++)
					{
						InvChanges[model.routes.get(temp.fromRoute).depot][p] -= temp.node.demandPerProduct[p];
						InvChanges[model.routes.get(temp.toRoute).depot][p]   += temp.node.demandPerProduct[p];
					}
					temp = temp.predecessor;
				}
				for (int d=0; d < model.depots; d++)
					if ( !model.hasSufficientInventory(d, InvChanges[d]) )
						sufficientInv = false;
			}
			
			//1.) the destination route fulfills all constraints
			if (( bufferTour>=0 && bufferCapacity >= 0 && sufficientInv))
			{						
				//1. a) start a new chain, if there is no predecessor
				if ( bestRelocationOption == null )
				{					
					ArrayList<Integer> forbidden = new ArrayList<Integer>();
					ArrayList<Integer> destinations = new ArrayList<Integer>();
					EjectionChain.Relocation opt = new Relocation( relocationSaving_p, relocationSaving, relocate, oldRoute, newRoute, insertionPos, 0, bufferCapacity, bufferTour, oldDetour, newDetour, null, forbidden, destinations  );
					opt.relocateDestinations.add(1000*newRoute+insertionPos);
					opt.relocateDestinations.add(1000*oldRoute+relocate.routePosition);
					opt.relocateDestinations.add(1000*oldRoute+relocate.routePosition+1);
					opt.forbiddenNodes.add(relocate.index);
					opt.forbiddenNodes.add(relocate.neighbours[0].index);
					opt.forbiddenNodes.add(relocate.neighbours[1].index);
					
					EjectionChain.Relocation valid = opt.copy();
					valid.forbiddenNodes.add(model.routes.get(newRoute).nodes.get(insertionPos).index);
					valid.forbiddenNodes.add(model.routes.get(newRoute).nodes.get(insertionPos-1).index);		
					addvalidChain( valid);

					relocations_all.add(opt);
					relocations_to_be_extended.add(relocations_all.size()-1);
				}
				//1. b) continue an existing chain
				else
				{
					Relocation opt = new Relocation( relocationSaving_p, relocationSaving, relocate, oldRoute, newRoute, insertionPos, bestRelocationOption.chainPos + 1, bufferCapacity, bufferTour, oldDetour, newDetour, bestRelocationOption, bestRelocationOption.forbiddenNodes, bestRelocationOption.relocateDestinations);
					opt.forbiddenNodes.add(relocate.index);
					opt.forbiddenNodes.add(relocate.neighbours[0].index);
					opt.forbiddenNodes.add(relocate.neighbours[1].index);
					opt.forbiddenNodes.add(model.routes.get(bestRelocationOption.toRoute).nodes.get(bestRelocationOption.newpos).index);
					opt.forbiddenNodes.add(model.routes.get(bestRelocationOption.toRoute).nodes.get(bestRelocationOption.newpos-1).index);
					opt.relocateDestinations.add(1000*newRoute+insertionPos);
					opt.relocateDestinations.add(1000*oldRoute+relocate.routePosition);
					opt.relocateDestinations.add(1000*oldRoute+relocate.routePosition+1);
					
					EjectionChain.Relocation valid = opt.copy();
					valid.forbiddenNodes.add(model.routes.get(newRoute).nodes.get(insertionPos).index);
					valid.forbiddenNodes.add(model.routes.get(newRoute).nodes.get(insertionPos-1).index);		
					addvalidChain( valid);
					
					relocations_all.add(opt);
					relocations_to_be_extended.add(relocations_all.size()-1);
				}							
			}//buffer ok
			//2.) the destination route is not feasible
			else
			if (true)//lastlevel
			{
				//2. a) start a new chain, if there is no predecessor
				if ( bestRelocationOption == null)
				{
					ArrayList<Relocation> opt = new ArrayList<Relocation>();
					ArrayList<Integer> forbidden = new ArrayList<Integer>();
					ArrayList<Integer> destinations = new ArrayList<Integer>();				
					Relocation relOpt =	new Relocation( relocationSaving_p, relocationSaving, relocate, oldRoute, newRoute, insertionPos, 0, bufferCapacity, bufferTour, oldDetour, newDetour, null, forbidden, destinations );
					relOpt.relocateDestinations.add(1000*oldRoute+relocate.routePosition);
					relOpt.relocateDestinations.add(1000*newRoute+insertionPos);
					relOpt.relocateDestinations.add(1000*oldRoute+relocate.routePosition+1);
					relOpt.forbiddenNodes.add(relocate.index);
					relOpt.forbiddenNodes.add(relocate.neighbours[0].index);
					relOpt.forbiddenNodes.add(relocate.neighbours[1].index);
					opt.add( relOpt );
					
					relocations_all.add(relOpt);
					relocations_to_be_extended.add(relocations_all.size()-1);					
				}
				//2. b) continue an existing chain
				else
				{
					Relocation opt = new Relocation( relocationSaving_p, relocationSaving, relocate, oldRoute, newRoute, insertionPos, bestRelocationOption.chainPos + 1, bufferCapacity, bufferTour, oldDetour, newDetour, bestRelocationOption, bestRelocationOption.forbiddenNodes, bestRelocationOption.relocateDestinations);
					opt.relocateDestinations.add(1000*newRoute+insertionPos);
					opt.relocateDestinations.add(1000*oldRoute+relocate.routePosition);
					opt.relocateDestinations.add(1000*oldRoute+relocate.routePosition+1);
					opt.forbiddenNodes.add(relocate.index);
					opt.forbiddenNodes.add(relocate.neighbours[0].index);
					opt.forbiddenNodes.add(relocate.neighbours[1].index);
					opt.forbiddenNodes.add(model.routes.get(bestRelocationOption.toRoute).nodes.get(bestRelocationOption.newpos).index);
					opt.forbiddenNodes.add(model.routes.get(bestRelocationOption.toRoute).nodes.get(bestRelocationOption.newpos-1).index);
					
					relocations_all.add(opt);
					relocations_to_be_extended.add(relocations_all.size()-1);					
				}
			}
	}//function
	
	
	private void addvalidChain(Relocation opt)
	{
		ArrayList<Relocation> validChain = new ArrayList<Relocation>();
		Relocation next = opt;
		while (next != null)
		{
			validChain.add(next.copy());
			next = next.predecessor;
		}
		candidate_moves.add(validChain);
	}
	
	//for all nodes to start from
	//is original route feasible if removed?
	//get nearest neighbours
	//is overall improvement?

}//CLASS
