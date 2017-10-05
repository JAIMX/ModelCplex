package version6;

import java.io.IOException;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;
import version6.Data.Edge;

public class Bender {

	/**
	 * Minimize Z=c'*y+f'x s.t A*y+B*x>=b y>=0 X in X
	 */

	static final double FUZZ = 1.0e-7;

	IloCplex master;
	IloCplex sub;
	IloNumVar[] y;
//	IloIntVar[] x;
//	IloNumVar z;
	double[] c, f, b;
	double[][] A, B;
	int numX, numY, numConstraint;
	double[] xValues, yValues; // subproblem y values ||master x values
	IloRange[] subConstraint;
	// HashMap<IloConstraint, Integer> rhs; // here we use rhs to record the index
	// of constraints in sub
	Data data;
	double tolerance, UB, LB, zMaster;
	
	//for master
	ArrayList<double[]> feasibleCut, optimalCut;
	double[] optPrice;
	double[] feaPrice;
	double[] truckPrice;
	double[] initialX;

	public Bender(double[] c, double[] f, double[] b, double[][] A, double[][] B, Data data, double tolerance,double[] initialX)
			throws IloException {
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
		
		this.initialX=initialX;

		// rhs = new HashMap<IloConstraint, Integer>();

		// set up the master problem(which initially has no constraint)
		master = new IloCplex();
		BuildMaster();
		master.setOut(null);

		
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
//		x = new IloIntVar[numX];
//		for (int k = 0; k < data.numberOfTrucks; k++) {
//			for (int edgeIndex = 0; edgeIndex < data.edgeSet.size(); edgeIndex++) {
//				Edge tempedge = data.edgeSet.get(edgeIndex);
//				x[k * data.edgeSet.size() + edgeIndex] = master.intVar(0, 1,
//						"x" + tempedge.u + "," + tempedge.t1 + "," + tempedge.v + "," + tempedge.t2 + "," + k);
//			}
//		}
//
//		// x = master.intVarArray(numX, 0, 1);
//		z = master.numVar(Double.MIN_VALUE, Double.MAX_VALUE, "ArificialVar");
//		master.addMinimize(master.sum(z, master.scalProd(f, x)), "Obj");
		
		//record the parameter of constraints in master
		feasibleCut = new ArrayList<double[]>();
		optimalCut = new ArrayList<double[]>();
		

	}

	public final void solve() throws IloException {
//		master.solve();
//		zMaster = master.getValue(z);
//		xValues = master.getValues(x);
		
		xValues=initialX;

		double tempConst = 0;
		for (int i = 0; i < numX; i++) {
			tempConst += f[i] * xValues[i];
		}

		UB = Double.MAX_VALUE;
		LB = Double.MIN_VALUE;

		while (UB - LB > tolerance) {
			
			System.out.println("Now the upper bound= " + UB);
			System.out.println("lower bound= " + LB);
			// System.out.print("currentX= ");
			// System.out.println(Arrays.toString(xValues));
			System.out.println();

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
//			IloNumExpr expr = master.numExpr();

			// System.out.println("The subproblem is "+status.toString());

			if (status == IloCplex.Status.Infeasible) {
				// subproblem is infeasible -- add a feasibility cut
				// first step: get a Farkas certificate, corresponding to a dual ray
				// along which the dual is unbounded

				IloConstraint[] constraints = new IloConstraint[numConstraint];
				double[] coefficients = new double[numConstraint];
				sub.dualFarkas(constraints, coefficients);

				double[] mu = coefficients;
				double[] tempFeaCut=new double[numX+1];

				// add a feasibility cut
				double temp = 0;
				for (int i = 0; i < numConstraint; i++) {
					temp += mu[i] * b[i];
				}
				tempFeaCut[numX]=temp;

				for (int i = 0; i < numX; i++) {
					double para = 0;
					for (int j = 0; j < numConstraint; j++) {
						para += mu[j] * B[j][i];
					}
					para = -para;
					
					tempFeaCut[i]=para;
//					expr = master.sum(expr, master.prod(para, x[i]));
				}

				feasibleCut.add(tempFeaCut);
//				IloConstraint r = master.addGe(0, expr);
				// System.out.println("\n>>> Adding feasibility cut: " + r + "\n");
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
				double[] tempOptCut=new double[numX+1];
				
				double temp = 0;
				for (int i = 0; i < numConstraint; i++) {
					temp += mu[i] * b[i];
				}
				tempOptCut[numX]=temp;

				for (int i = 0; i < numX; i++) {
					double para = 0;
					for (int j = 0; j < numConstraint; j++) {
						para += mu[j] * B[j][i];
					}
					para = -para;
					
					tempOptCut[i]=para;
//					expr = master.sum(expr, master.prod(para, x[i]));
				}

				optimalCut.add(tempOptCut);

//				IloConstraint r = master.addGe(z, expr);
				// System.out.println("\n>>> Adding optimality cut: " + r + "\n");
				System.out.println("\n>>> Adding optimality cut: " + "\n");

			} else {
				// unexpected status -- report but do nothing
				System.err.println("\n!!! Unexpected subproblem solution status: " + status + "\n");
			}

			xValues=BMP();
//			xValues = master.getValues(x);
			tempConst = 0;
			for (int i = 0; i < numX; i++) {
				tempConst += f[i] * xValues[i];
			}
			// LB = master.getValue(z)+tempConst;
			LB = master.getObjValue();

		}
		
		System.out.println("Now the upper bound= " + UB);
		System.out.println("lower bound= " + LB);
		System.out.println();
		
		
		
		
		
		if (master.solve()) {
			System.out.println("optimal obj= " + master.getObjValue());
			// double[] xValues = master.getValues(x);
			// System.out.println("x= " + Arrays.toString(xValues));
			// System.out.println("y= " + Arrays.toString(yValues));
		} else {
			System.out.println("The last master's status is " + master.getStatus().toString());
		}
		
		
		
		
		
		
	}

	public static void main(String[] args) throws IOException, IloException {
		Data data = new Data();
//		 data.readData("./data/temp.txt");
		// System.out.println("Read data done!");
//		 data.readData("./data/out_small.txt");
//		data.readData("./data/data1.txt");
		data.readData("./data/data2.txt");
		data.graphTransfer();
		// System.out.println("Graph transfer done!");
		data.matrixGenerator();
		// System.out.println("MatrixGenerator done!");
		double tolerance = 0;

		Bender test = new Bender(data.c, data.f, data.bb, data.A, data.B, data, tolerance,data.generateInitialx());
	}

}
