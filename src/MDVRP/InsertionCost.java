package MDVRP;

public class InsertionCost implements Comparable<InsertionCost>
{
		public double costs_p;
		public double costs;
		public int route;
		public int position;
			
		public InsertionCost(double c1, double c2, int r, int p)
		{
			costs_p = c1;
			costs 	= c2;
			route	= r;
			position= p;
		}
		
		public InsertionCost copy()
		{
			InsertionCost copy = new InsertionCost(costs_p, costs, route, position );
			return copy;
		}
		
		//TODO copy


		@Override
		public int compareTo(InsertionCost i) 
		{
			if(i.costs_p<this.costs_p)
				return 1;
			else if(i.costs_p == this.costs_p)
				return 0;
			else
				return -1;
		}
}
