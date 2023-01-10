package MDVRP;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class ClarkWright {
	
	private  class Saving implements Comparable<Saving>
	{
		public double saving;
		public int from;
		public int to;
			
		public Saving(double s, int f, int t)
		{
			saving = s;
			from = f;
			to = t;
		}

		@Override
		public int compareTo(Saving o) 
		{
			if(o.saving<this.saving)
				return -1;
			else if(o.saving == this.saving)
				return 0;
			else
				return 1;
		}

	}//class

	@SuppressWarnings("rawtypes")
	 void clarkWrightParallel( MDVRPModel model, int depot )
	{
		boolean[] visited	= new boolean[model.customers];
		int depot_node	= model.customers + depot;
		
		//---------------------------------------------
		//compute and order the savings for all pairs of customers
		//---------------------------------------------
		ArrayList<Saving> savings = new ArrayList<Saving>();
		for (int i=0; i<model.customers; i++)
		if (model.nodes[i].depot == depot)
		{
			Set<Entry<Integer, Float>> nearestDistances = model.nodes[i].distances_for_CW.entrySet();
			Iterator<Entry<Integer, Float>> iterator 	= nearestDistances.iterator();
		    while(iterator.hasNext()) 
			//for (int j=i+1; j<model.customers; j++)
		    {
		    	Map.Entry mentry = (Map.Entry)iterator.next();
		    	int neighbour	= (int) mentry.getKey();
		    	//is the neighbour a customer that is allocated to the same depot 
				if (model.nodes[neighbour].depot == depot && neighbour < model.customers )
				{
					double saving	 = model.getDistance(i, depot_node) + model.getDistance(neighbour, depot_node) - model.getDistance(i, neighbour);
					if (saving < 0.1)
						saving = 0.1;
					Saving s = new Saving( saving, i, neighbour );
					savings.add(s);
				}
			}
		    //System.out.print(i+" ");
		}
		else
			visited[i] = true;
		Collections.sort( savings );
		 //System.out.println();
				
		//---------------------------------------------
		//iteratively add edges between pairs of nodes with the highest savings
		//---------------------------------------------		
		int countVisits		= 0;
		//int countNewRoutes	= 0;
		boolean[] extensionPoints	= new boolean[model.customers];
		boolean[] interiorPoints	= new boolean[model.customers];		
		int index_first_route		= model.routes.size();
		int index_current_route 	= index_first_route;
		Iterator<Saving> itr = savings.iterator();
		
		//2E LRP
		/*boolean splitDeliveryPossible = false;
		if ( model.isFirstLevel )
		{
			for (Node n: model.nodes)
				if ( n.demand >= model.capacityLimitVehicle )
				{
					Route r = new Route();
					r.nodes.add(model.nodes[depot_node]);
					r.nodes.add(n);
					r.nodes.add(model.nodes[depot_node]);
					r.length 			= 1;
					r.load				= model.capacityLimitVehicle;
					r.depot				= depot;
					model.routes.add( r.copy() );
					n.demand -= model.capacityLimitVehicle;
				}
			//is a split delivery possible?
			double sumInventory = 0;
			for (Node n: model.nodes)
				sumInventory += n.demand;
			
			if ( Math.ceil( (double) sumInventory / (double) model.capacityLimitVehicle ) < model.customers )
				splitDeliveryPossible = true;
			splitDeliveryPossible = false;
			//for (Node n: model.nodes)
			//System.out.print(n.demand + " ");
			//System.out.println( sumInventory + " " + model.capacityLimitVehicle);
		}
		ArrayList<Integer> splitDeliveryNodes = new ArrayList<Integer>();
		boolean splittedBefore = false;//split only once
		*/
		
			
		while (itr.hasNext() )
		{	
			//---------------------------------------------
			// find the next edge
			//---------------------------------------------
			Saving s = itr.next();
				int n1 = s.from;
				int n2 = s.to;
				//System.out.println(n1 + " " + n2 + " " + s.saving)
				//---------------------------------------------
				// case 1: one of the nodes is an interior Point
				//---------------------------------------------
				if (interiorPoints[n1] || interiorPoints[n2])
				{
				}
				//---------------------------------------------
				// case 2: both nodes are neither extension nor interior points
				//---------------------------------------------
				else if (!extensionPoints[n1] && !extensionPoints[n2] )
				{
					if ( model.nodes[ n1 ].demand + model.nodes[ n2 ].demand <= model.capacityLimitVehicle )//split delivery not allowed?
					{
						//start a new route with these two nodes
						Route r = new Route();
						r.nodes.add(model.nodes[depot_node]);
						r.nodes.add(model.nodes[ n1 ]);
						r.nodes.add(model.nodes[ n2 ]);
						r.nodes.add(model.nodes[depot_node]);
						r.length 			= 2;
						r.load				= model.nodes[ n1 ].demand + model.nodes[ n2 ].demand;
						r.depot				= depot;
						r.cost				= model.getDistance(n1, n2) + model.getDistance(n1, depot_node) + model.getDistance(n2, depot_node);//model.nodes[ n1 ].distances[ n2 ] + model.nodes[ n1 ].distances[ depot ] + model.nodes[ n2 ].distances[ depot ];
						countVisits = countVisits + 2;
						extensionPoints[n1]	=true;
						extensionPoints[n2]	=true;
						visited[n1]=true;
						visited[n2]=true;
						model.nodes[ n1 ].routeIndex = index_current_route;
						model.nodes[ n2 ].routeIndex = index_current_route;						
						model.routes.add( r.copy() );
						index_current_route++;	
						//2E LRP
						/*if ( model.nodes[ n1 ].demand + model.nodes[ n2 ].demand > model.capacityLimitVehicle || splitDeliveryNodes.contains(n1) || splitDeliveryNodes.contains(n2) )
						{
							extensionPoints[n1]	=false;
							extensionPoints[n2]	=false;
							splitDeliveryNodes.add(n1);
							splitDeliveryNodes.add(n2);
							if ( splittedBefore )
								for (Integer n : splitDeliveryNodes)
									extensionPoints[n] = true;
							else
								splittedBefore = true;						
						}*/
							
					}
				}
				//---------------------------------------------
				// case 3a: the first node is an extension point 
				//---------------------------------------------
				else if (extensionPoints[n1] && !extensionPoints[n2])
				{
					Route r = model.routes.get( model.nodes[n1].routeIndex );
					double costs_add = model.getDistance(n1, n2) + model.getDistance(n2, depot_node) - model.getDistance(n1, depot_node);
					if ( r.load + model.nodes[ n2 ].demand <= model.capacityLimitVehicle && r.cost + costs_add <= model.maxTourLength)
					{						
						//extend the route
						if ( r.nodes.get(1).index == n1 )
							r.nodes.add(1, model.nodes[ n2 ]);
						else
							r.nodes.add(r.length+1, model.nodes[ n2 ]);
						r.length++;
						r.load	+= model.nodes[ n2 ].demand;
						r.cost += costs_add;
						model.nodes[ n2 ].routeIndex = model.nodes[n1].routeIndex;
						countVisits++;
						extensionPoints[n1] = false;
						extensionPoints[n2]	= true;
						interiorPoints[n1]	= true;
						visited[n2]=true;
					}
				}
				//---------------------------------------------
				// case 3b: the second node is an extension point 
				//---------------------------------------------
				else if (!extensionPoints[n1] && extensionPoints[n2])
				{
					Route r = model.routes.get( model.nodes[n2].routeIndex );
					double costs_added = model.getDistance(n2, n1) + model.getDistance(n1, depot_node) - model.getDistance(n2, depot_node);
					if ( r.load + model.nodes[ n1 ].demand <= model.capacityLimitVehicle && r.cost + costs_added <= model.maxTourLength)
					{						
						//extend the route
						if ( r.nodes.get(1).index == n2 )
							r.nodes.add(1, model.nodes[ n1 ]);
						else
							r.nodes.add(r.length+1, model.nodes[ n1 ]);
						r.length++;
						r.load	+= model.nodes[ n1 ].demand;
						r.cost += costs_added;
						model.nodes[ n1 ].routeIndex = model.nodes[n2].routeIndex;
						countVisits++;
						extensionPoints[n2]	= false;
						extensionPoints[n1]	= true;
						interiorPoints[n2]	= true;
						visited[n1]=true;	
					}				
				}
				//---------------------------------------------
				// case 4: both nodes are extension points; merge the two routes
				//---------------------------------------------
				else if (extensionPoints[n1] && extensionPoints[n2])
				{
					if (model.nodes[n1].routeIndex != model.nodes[n2].routeIndex )
					{
						Route r1 = model.routes.get( model.nodes[n1].routeIndex );
						Route r2 = model.routes.get( model.nodes[n2].routeIndex );
						double costs_merge = r1.cost + r2.cost + model.getDistance(n2, n1) - model.getDistance(n1, depot_node) - model.getDistance(n2, depot_node);
						if ( r1.load + r2.load <= model.capacityLimitVehicle && costs_merge <= model.maxTourLength )
						{
							//merge the routes and safe as new route
							//bring into format d -...- n1 - n2 -...- d
							if ( r1.nodes.get(1).index == n1 )
								Collections.reverse(r1.nodes);
							if ( r2.nodes.get(1).index != n2 )
								Collections.reverse(r2.nodes);
							//remove two depots
							r1.nodes.remove( r1.nodes.size()-1 );
							r2.nodes.remove( 0 );
							//attach r2 to r1
							r1.nodes.addAll(r1.nodes.size(), r2.nodes);
							r1.length 	+= r2.length;
							r1.load		+= r2.load;
							Route temp	= r1.copy();
							
							//update the routes
							//adapt the routeIndices of all other routes
							int index1 = model.nodes[n1].routeIndex;
							int index2 = model.nodes[n2].routeIndex;
							for (int i=index1+1; i<model.routes.size(); i++)
								for (int j=1; j<model.routes.get(i).length+1; j++)
									model.routes.get(i).nodes.get(j).routeIndex--;
							for (int i=index2+1; i<model.routes.size(); i++)
								for (int j=1; j<model.routes.get(i).length+1; j++)
									model.routes.get(i).nodes.get(j).routeIndex--;
							//remove the two previous routes
							if ( index1 < index2 )
							{
								model.routes.remove(index2);
								model.routes.remove(index1);
							}
							else
							{
								model.routes.remove(index1);
								model.routes.remove(index2);
							}	
							//insert the merged route
							for (int i=1; i<temp.length+1; i++)
								temp.nodes.get(i).routeIndex = index_current_route-2;
							temp.cost = costs_merge;
							model.routes.add(temp.copy());
							index_current_route--;
	
							extensionPoints[n1]	= false;
							extensionPoints[n2]	= false;
							interiorPoints[n1]	= true;
							interiorPoints[n2]	= true;
						}
					}
				}
		}//while not all visited
		//catch the case that not all customers fitted in the routes
		if ( countVisits < model.customers)
		{
			for (int i=0; i<model.customers; i++)
				if (!visited[i] && model.nodes[i].depot == depot )
				{
					//create an own route for this customer
					Route r = new Route();
					r.nodes.add(model.nodes[depot_node]);
					r.nodes.add(model.nodes[ i ]);
					r.nodes.add(model.nodes[depot_node]);
					r.length 			= 1;
					r.load				= model.nodes[ i ].demand;
					r.depot				= depot;
					model.nodes[ i ].routeIndex = index_current_route;
					countVisits++;
					model.routes.add(r.copy());
				}
		}
		//set route and position values, and update cost values and depots;
		for (int i=index_first_route; i<model.routes.size(); i++)
		{
			model.routes.get(i).cost = model.getDistance(model.routes.get(i).nodes.get(0).index, model.routes.get(i).nodes.get(1).index);
			for (int j=1; j<model.routes.get(i).length+1; j++)
			{
				model.routes.get(i).nodes.get(j).routeIndex = i;
				//System.out.print(model.routes.get(i).nodes.get(j).index + ",");
				model.routes.get(i).nodes.get(j).routePosition = j;
				model.routes.get(i).cost += model.getDistance(model.routes.get(i).nodes.get(j).index, model.routes.get(i).nodes.get(j+1).index);//model.routes.get(i).nodes.get(j).distances[model.routes.get(i).nodes.get(j+1).index];
			}
			//System.out.println();
		}
		//create empty routes (which might be filled during the heuristic)
		for (int i=0; i<1; i++)
		{
			Route r = new Route();
			r.nodes.add(model.nodes[depot_node]);
			r.nodes.add(model.nodes[depot_node]);
			r.length 			= 0;
			r.load				= 0;
			r.depot				= depot;
			model.routes.add(r.copy());
		}
		//System.out.println("CW " + model.computeSolutionCosts());
	}
	
	 public void clarkWright_MDVRP( MDVRPModel model )
	{
		//each customer is allocated to its closest depot, and then routes are computed with CW per depot
		for (int cust=0; cust<model.customers; cust++)
		{
			double shorted_dist = 999999;
			for (int depot=0; depot<model.depots; depot++)
			{				
				double dist = model.getDistance(cust, model.customers + depot);
				if ( dist < shorted_dist)
				{
					model.nodes[cust].depot = depot;
					shorted_dist = dist;
				}
			}
		}		
		//construct initial routes, for each depot separately
		for (int d=0; d < model.depots; d++)
			clarkWrightParallel( model, d );
	}
	
	 public boolean clarkWright_MPMDVRPI( MDVRPModel model, ArrayList<Integer> pool, int poolLimit, boolean suppressInfeasibility )
	{
		//if (!pool.isEmpty() && !suppressInfeasibility)
		//	;//System.out.println(poolLimit + " - " +  (int)(model.depotInventory[pool.get(0)][0]+model.depotInventory[pool.get(1)][0]) );
		//pool.clear();
		
		//-------------------------------------------------------------------------
		//regret heuristic - for each customer, the distance between to its nearest depot and second nearest depot is computed.
		//the difference between both distances is saved as regret.
		//iteratively, the customers with the highest regret are allocated to its nearest depot
		//if that still has sufficient inventory
		//-------------------------------------------------------------------------
		//reset the inventory of all depots
		model.remainingInventory = new int[model.depots][model.products];
		for (int d=0; d<model.depots; d++)
			for (int p=0; p<model.products; p++)
				model.remainingInventory[d][p] = model.depotInventory[d][p];
				
		//compute regrets for each customer
		ArrayList<Saving> regret_list = new ArrayList<Saving>();
		for (int cust=0; cust<model.customers; cust++)
		{
			double dist_nearest_depot = 999999;
			double dist_second_nearest_depot = 999999;
			int nearestDepot = -1;
			Node node = model.nodes[cust];
			for (int depot=0; depot<model.depots; depot++)
			{
				//can the customer be delivered by that depot?
				if ( model.hasSufficientInventory(depot, node.demandPerProduct) )
				{
					double dist = model.getDistance(cust, model.customers+depot);
					if ( dist < dist_nearest_depot)
					{
						dist_second_nearest_depot 	= dist_nearest_depot;
						dist_nearest_depot 			= dist;
						nearestDepot 				= depot;
					}
					else if (dist < dist_second_nearest_depot)
					{
						dist_second_nearest_depot = dist;
					}
					
				}
			}
			if ( nearestDepot == -1 )
				System.err.println("Inventory constraints too tight. One customer cannot be supplied by any depot.");
			Saving s = new Saving( dist_second_nearest_depot - dist_nearest_depot, nearestDepot, node.index );
			regret_list.add(s);
		}
		Collections.sort( regret_list );
		//on the basis of regrets, assign customers to depots

		if (!assignCustomersToDepot( regret_list, model ) && false )
			return false;
		
		//assignCustomersToDepotSimple( regret_list, model );

		for (int d=0; d < model.depots; d++)
			clarkWrightParallel( model, d );
		return true;
		//compute routes with CW
		
	}
	
	
	 private boolean assignCustomersToDepot(ArrayList<Saving> regret_list, MDVRPModel model)
	{
		boolean feasibleAlloction = true;
		//each customer has no assigned depot
		for (int i=0; i<model.customers; i++)
			model.nodes[i].depot = -1;
		
		while ( !regret_list.isEmpty() )
		{
			//assign the customer with the highest regret
			Saving s = regret_list.get(0);
			int assignedDepot = s.from;
			model.nodes[s.to].depot = assignedDepot;

			for (int p=0; p<model.products; p++)
				model.remainingInventory[assignedDepot][p] -= model.nodes[s.to].demandPerProduct[p];
			regret_list.remove(0);
			
			
			//update the regrets of all unassigned customers			
			for (int i=0; i<regret_list.size(); i++)
			{
				int customer = regret_list.get(i).to;
				//recompute regrets, if the primary depot no longer has sufficient inventory
				if ( regret_list.get(i).from == assignedDepot )
					if ( !model.hasSufficientInventory(assignedDepot, model.nodes[customer].demandPerProduct) )
					{
						double dist_nearest_depot 		= 999999;
						double dist_second_nearest_depot= 999999;
						int nearestDepot 				= -1;						
						for (int d=0; d<model.depots; d++)
							if (model.hasSufficientInventory(d, model.nodes[customer].demandPerProduct) ){
								double dist = model.getDistance(customer, model.customers+d);
								if ( dist < dist_nearest_depot)
								{
									dist_second_nearest_depot = dist_nearest_depot;
									dist_nearest_depot = dist;
									nearestDepot = d;
								}
								else if (dist < dist_second_nearest_depot)
								{
									dist_second_nearest_depot = dist;
								}	
							}
						
						if ( nearestDepot == -1 )
						{
							//System.err.println("WARNING: Cannot obtain initial feasible solution with the available inventory levels.");
							//model.nodes[customer].depot = 0;
							//model.remainingInventory[0][0]	-= model.nodes[customer].demand;
							//regret_list.remove(i);
							//i--;
							regret_list.get(i).saving = 0;
							feasibleAlloction = false;
						}
						else	
						{
							Saving replace = new Saving( dist_second_nearest_depot-dist_nearest_depot, nearestDepot, customer );
							regret_list.remove(i);
							regret_list.add(i,replace);
						}
					}
			}
			Collections.sort( regret_list );			
		}//allocate all customers
		
		return feasibleAlloction;
	}
	
	
	 private boolean assignCustomersToDepotSimple(ArrayList<Saving> regret_list, MDVRPModel model)
	{
		int inventoryNeededInPool = 0;
		boolean feasibleAlloction = true;
		//each customer has no assigned depot
		for (int i=0; i<model.customers; i++)
			model.nodes[i].depot = -1;
		
		while ( !regret_list.isEmpty() )
		{
			//assign the customer with the highest regret
			Saving s = regret_list.get(0);
			int assignedDepot = s.from;
			if (model.hasSufficientInventory(assignedDepot, model.nodes[s.to].demandPerProduct) )
					;
			else 
				for (int d=0; d<model.depots; d++)
					if (model.hasSufficientInventory(d, model.nodes[s.to].demandPerProduct) )
						assignedDepot = d;

			model.nodes[s.to].depot = assignedDepot;
			for (int p=0; p<model.products; p++)
				model.remainingInventory[assignedDepot][p] -= model.nodes[s.to].demandPerProduct[p];
			regret_list.remove(0);			
				
		}//allocate all customers
		
		return feasibleAlloction;
	}
	
}
