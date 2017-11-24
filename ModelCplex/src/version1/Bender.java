package version1;

import java.io.IOException;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;
import version1.Data.Edge;

public class Bender {

	/**
	 * Minimize Z=c'*y+f'x s.t A*y+B*x>=b y>=0 X in X
	 */

	static final double FUZZ = 1.0e-7;

	IloCplex master;
	IloCplex sub;
	IloNumVar[] y;
	IloIntVar[] x;
	IloNumVar z;
	double[] c, f, b;
	double[][] A, B;
	int numX, numY, numConstraint;
	double[] xValues, yValues; // subproblem y values ||master x values
	IloRange[] subConstraint;
	// HashMap<IloConstraint, Integer> rhs; // here we use rhs to record the index
	// of constraints in sub
	Data data;
	double tolerance, UB, LB, zMaster;
	
	double startTime, currentTime, lastTime;
	double time1;//subproblem total time
	double time2;//extreme ray or extreme point total time
	double time3;//master problem total time
	

	public Bender(double[] c, double[] f, double[] b, double[][] A, double[][] B, Data data, double tolerance)
			throws IloException {
		startTime = System.currentTimeMillis();
		lastTime=startTime;
		time1=0;
		time2=0;
		time3=0;
		// System.out.println("c= "+Arrays.toString(c));
		// System.out.println("f= "+Arrays.toString(f));
		// System.out.println("b= "+Arrays.toString(b));
		// System.out.println("f= "+Arrays.toString(f));

		numY = c.length;
		numX = f.length;
		numConstraint = A.length;
		this.data = data;
		this.tolerance = tolerance;

		// this.c = Arrays.copyOf(c, c.length);
		// this.f = Arrays.copyOf(f, f.length);
		// this.b = Arrays.copyOf(b, b.length);
		this.c = c;
		this.f = f;
		this.b = b;

		this.A = A;
		this.B = B;

		// rhs = new HashMap<IloConstraint, Integer>();

		// set up the master problem(which initially has no constraint)
		master = new IloCplex();
		BuildMaster();
		master.setOut(null);
		// System.out.println("Initial master is "+master.toString());
		// master.use(new BendersCallback());

		// set up the subproblem
		sub = new IloCplex();
		BuildSub();

		solve();

	}

	public void BuildSub() throws IloException {

		y = new IloNumVar[numY];
		for (int p = 0; p < data.numberOfDemandPair; p++) {
			for (int edgeIndex = 0; edgeIndex < data.numOfEdge12; edgeIndex++) {
				Edge tempedge = data.edgeSet.get(edgeIndex);
				y[p * data.numOfEdge12 + edgeIndex] = sub.numVar(0, Double.MAX_VALUE,
						"y" + tempedge.u + "," + tempedge.t1 + "," + tempedge.v + "," + tempedge.t2 + "," + p);
			}
		}
		// y = sub.numVarArray(numY, 0, Double.MAX_VALUE);

		sub.addMinimize(sub.scalProd(c, y));
		// constrain to be satisfied -- record the constraints for use later
		subConstraint = new IloRange[numConstraint];
		IloLinearNumExpr expr = sub.linearNumExpr();

		// add constraints (initially all zero, which makes the subproblem
		// infeasible)
		// -- record the constraints for use later
		// -- also map each constraint to the corresponding binary variable in the
		// master (for decoding Farkas certificates)
		for (int row = 0; row < numConstraint; row++) {
			expr.clear();
			expr = sub.scalProd(A[row], y);
			subConstraint[row] = sub.addGe(expr, 0, "subConstraint_" + row);
			// rhs.put(subConstraint[row], master.diff(b[row], master.scalProd(B[row], x)));
			// rhs.put(subConstraint[row], row);
		}

		// System.out.println("initial subproblem is "+sub.toString());
		// disable presolving of the subproblem (if the presolver realizes the
		// subproblem is infeasible, we do not get a dual ray)
		sub.setParam(IloCplex.BooleanParam.PreInd, false);
		// force use of the dual simplex algorithm to get a Farkas certificate
		sub.setParam(IloCplex.IntParam.RootAlg, 2);
		// suppress subproblem output to reduce clutter
		sub.setOut(null);
	}

	public void BuildMaster() throws IloException {
		x = new IloIntVar[numX];
		for (int k = 0; k < data.numberOfTrucks; k++) {
			for (int edgeIndex = 0; edgeIndex < data.edgeSet.size(); edgeIndex++) {
				Edge tempedge = data.edgeSet.get(edgeIndex);
				x[k * data.edgeSet.size() + edgeIndex] = master.intVar(0, 1,
						"x" + tempedge.u + "," + tempedge.t1 + "," + tempedge.v + "," + tempedge.t2 + "," + k);
			}
		}

		// x = master.intVarArray(numX, 0, 1);
		z = master.numVar(Double.MIN_VALUE, Double.MAX_VALUE, "ArificialVar");
		master.addMinimize(master.sum(z, master.scalProd(f, x)), "Obj");
		// master.addMinimize(master.scalProd(f, x), "Obj");

		// constraints about road construction
		// ---constraint 1-4---//
		for (int k = 0; k < data.numberOfTrucks; k++) {

			// ---constraint 1---//
			for (int i = 0; i < data.numberOfCities * (data.T + 1); i++) {
				IloLinearNumExpr constraint1 = master.linearNumExpr();
				for (int edgeIndex : data.distance.get(i)) {
					constraint1.addTerm(1, x[k * data.edgeSet.size() + edgeIndex]);
				}

				for (int edgeIndex : data.distanceReverse.get(i)) {
					constraint1.addTerm(-1, x[k * data.edgeSet.size() + edgeIndex]);
				}

				master.addEq(0, constraint1);
			}

			// ---constraint 2---//
			IloLinearNumExpr constraint2 = master.linearNumExpr();
			int ok = data.numberOfCities * (data.T + 1) + data.truckStartNode[k];
			for (int edgeIndex : data.distance.get(ok)) {
				constraint2.addTerm(1, x[k * data.edgeSet.size() + edgeIndex]);
			}
			master.addGe(1, constraint2);

			// ---constraint 3---//
			IloLinearNumExpr constraint3 = master.linearNumExpr();

			for (int ok2 = data.numberOfCities * (data.T + 1); ok2 < data.numberOfCities * (data.T + 2); ok2++) {
				if (ok2 != ok) {
					for (int edgeIndex : data.distance.get(ok2)) {
						constraint3.addTerm(1, x[k * data.edgeSet.size() + edgeIndex]);
					}
				}
			}
			master.addEq(0, constraint3);

			// ---constraint 4---//
			IloLinearNumExpr constraint4 = master.linearNumExpr();
			int dk = data.numberOfCities * (data.T + 2) + data.truckStartNode[k];

			for (int dk2 = data.numberOfCities * (data.T + 2); dk2 < data.numberOfCities * (data.T + 3); dk2++) {
				if (dk2 != dk) {
					for (int edgeIndex : data.distanceReverse.get(dk2)) {
						constraint4.addTerm(1, x[k * data.edgeSet.size() + edgeIndex]);
					}
				}
			}
			master.addEq(0, constraint4);

			// ---constraint 5---//
			IloLinearNumExpr constraint5 = master.linearNumExpr();
			for (int edgeIndex = 0; edgeIndex < data.edgeSet.size(); edgeIndex++) {
				Edge tempedge = data.edgeSet.get(edgeIndex);
				if (tempedge.setIndex == 1) {
					constraint5.addTerm(1, x[k * data.edgeSet.size() + edgeIndex]);
				}
			}
			master.addGe(data.legLimit, constraint5);

			// ---constraint 6---//
			IloLinearNumExpr constraint6 = master.linearNumExpr();
			for (int edgeIndex = 0; edgeIndex < data.edgeSet.size(); edgeIndex++) {
				constraint6.addTerm(data.edgeSet.get(edgeIndex).length, x[k * data.edgeSet.size() + edgeIndex]);
			}
			master.addGe(data.distanceLimit, constraint6);

		}
	}

	//name of X variables
	public String getName(int index) {
		String string = "";
		// k
		int k=index/data.edgeSet.size();
		int edgeIndex=index%data.edgeSet.size();
		Edge edge=data.edgeSet.get(edgeIndex);
		
		string="X"+edge.start+","+edge.end+","+k;

		return string;
	}
	
	public final void solve() throws IloException {
		master.solve();
		zMaster = master.getValue(z);
		xValues = master.getValues(x);

		double tempConst = 0;
		for (int i = 0; i < numX; i++) {
			tempConst += f[i] * xValues[i];
		}

		UB = Double.MAX_VALUE;
		LB = zMaster + tempConst;

		while (UB - LB > tolerance) {
			
			System.out.println("Now the upper bound= " + UB);
			System.out.println("lower bound= " + LB);
			// System.out.print("currentX= ");
			// System.out.println(Arrays.toString(xValues));
			System.out.println();
			
			double turnTime=0;
///----------------------------------------------------subProblem start----------------------------------------------///
			System.out.println();
			System.out.println("Now solve the subProblem");
			// set the supply constraint right-hand sides in the subproblem
			for (int row = 0; row < numConstraint; row++) {
				double temp = 0;
				for (int i = 0; i < numX; i++) {
					temp += B[row][i] * xValues[i];
				}
				subConstraint[row].setLB(b[row] - temp);
			}

			// solve the subproblem
//			System.out.println("||----Now solve the subProblem----||");
//			System.out.println();

			sub.solve();
			IloCplex.Status status = sub.getStatus();
			IloNumExpr expr = master.numExpr();
			
			currentTime = System.currentTimeMillis();
			double tempTime=(double)((currentTime-lastTime)/1000);
			time1+=tempTime;
			turnTime+=tempTime;
			System.out.println("subProblem solve end, we use "+tempTime+"s");
			lastTime=currentTime;

			///----------------------------------------------------subProblem end-----------------------------------------------///
			
			// System.out.println("The subproblem is "+status.toString());

			
			///----------------------------------------------------cut generation start----------------------------------------------///
			System.out.println();
			System.out.println("we start cut generation");
			if (status == IloCplex.Status.Infeasible) {
				// subproblem is infeasible -- add a feasibility cut
				// first step: get a Farkas certificate, corresponding to a dual ray
				// along which the dual is unbounded

				IloConstraint[] constraints = new IloConstraint[numConstraint];
				double[] coefficients = new double[numConstraint];
				sub.dualFarkas(constraints, coefficients);

				double mu[] = coefficients;

				// for (int row = 0; row < numConstraint; row++) {
				// mu[rhs.get(constraints[row])] = coefficients[row];
				// }

				// add a feasibility cut
				double tempconst = 0;
				for (int i = 0; i < numConstraint; i++) {
					tempconst += mu[i] * b[i];
				}
				expr = master.sum(tempconst, expr);

				for (int i = 0; i < numX; i++) {
					double para = 0;
					for (int j = 0; j < numConstraint; j++) {
						para += mu[j] * B[j][i];
					}
					para = -para;
					expr = master.sum(expr, master.prod(para, x[i]));
				}

				IloConstraint r = master.addGe(0, expr);
//				 System.out.println("\n>>> Adding feasibility cut: " + r + "\n");
				System.out.println("\n>>> Adding feasibility cut: " + "\n");
			} else if (status == IloCplex.Status.Optimal) {

				if (sub.getObjValue() + tempConst + FUZZ < UB) {
					UB = sub.getObjValue() + tempConst;
					yValues = sub.getValues(y);
					System.out.println("Find a  better solution, and updata UB!!!");
					System.out.println();
				}

				// add an optimality cut
				double[] mu = sub.getDuals(subConstraint);
//				System.out.println("mu= "+Arrays.toString(mu));
				double tempconst = 0;
				for (int i = 0; i < numConstraint; i++) {
					tempconst += mu[i] * b[i];
				}
				expr = master.sum(tempconst, expr);

				for (int i = 0; i < numX; i++) {
					double para = 0;
					for (int j = 0; j < numConstraint; j++) {
						para += mu[j] * B[j][i];
					}
					para = -para;
					expr = master.sum(expr, master.prod(para, x[i]));
				}

				IloConstraint r = master.addGe(z, expr);
//				 System.out.println("\n>>> Adding optimality cut: " + r + "\n");
				System.out.println("\n>>> Adding optimality cut: " + "\n");

			} else {
				// unexpected status -- report but do nothing
				System.err.println("\n!!! Unexpected subproblem solution status: " + status + "\n");
			}
			
			currentTime = System.currentTimeMillis();
			tempTime=(double)((currentTime-lastTime)/1000);
			time2+=tempTime;
			turnTime+=tempTime;
			System.out.println("cut generation end, we use "+tempTime+"s");
			lastTime=currentTime;
			
			///----------------------------------------------------cut generation end----------------------------------------------///

			
			///-------------------master solve start------------------///
			System.out.println();
			System.out.println("we start to solve master problem");
			
			master.solve();
			xValues = master.getValues(x);
//			for(int i=0;i<numX;i++) {
//				if(xValues[i]>0) {
//					
//					System.out.println(getName(i)+"="+xValues[i]);
//				}
//			}
			
			currentTime = System.currentTimeMillis();
			tempTime=(double)((currentTime-lastTime)/1000);
			time3+=tempTime;
			turnTime+=tempTime;
			System.out.println("master solve end, we use "+tempTime+"s");
			lastTime=currentTime;
			
			///-------------------master solve end------------------///
			tempConst = 0;
			for (int i = 0; i < numX; i++) {
				tempConst += f[i] * xValues[i];
			}
			// LB = master.getValue(z)+tempConst;
			LB = master.getObjValue();
			
			System.out.println("This turn we totally use "+turnTime+"s");
			System.out.println();
		}
		
		System.out.println("Now the upper bound= " + UB);
		System.out.println("lower bound= " + LB);
		System.out.println();
		
		
		
		
		
//		System.out.println("hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh");
		if (master.solve()) {
			System.out.println("optimal obj= " + master.getObjValue());
			// double[] xValues = master.getValues(x);
			// System.out.println("x= " + Arrays.toString(xValues));
			// System.out.println("y= " + Arrays.toString(yValues));
		} else {
			System.out.println("The last master's status is " + master.getStatus().toString());
		}
		
		
		System.out.println("subProblem total time= "+time1+"s");
		System.out.println("cut generation total time= "+time2+"s");
		System.out.println("master problem total time= "+time3+"s");
		
		
		
		
	}

	public static void main(String[] args) throws IOException, IloException {
		Data data = new Data();
//		 data.readData("./data/temp.txt");
//		 data.readData("./data/out_small.txt");
//		data.readData("./data/data1.txt");
//		data.readData("./data/data2.txt");
//		data.readData("./data/out_small2.txt");
//		data.readData("./data/out_small3.txt");
//		data.readData("./data/data1_1.txt");
//		data.readData("./data/out2.txt");
		data.readData("./data/report4_4.txt");
		data.graphTransfer();
		// System.out.println("Graph transfer done!");
		data.matrixGenerator();
		// System.out.println("MatrixGenerator done!");
		double tolerance = 0;

		Bender test = new Bender(data.c, data.f, data.bb, data.A, data.B, data, tolerance);
	}

}
