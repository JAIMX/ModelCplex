package version4;

import java.io.IOException;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;
import version4.Data.Edge;

public class Bender {

	static final double FUZZ = 1.0e-7;

	IloCplex master;
	IloCplex sub;
	IloNumVar[] y;
	IloIntVar[] x;
	IloNumVar z;
	double[] c, f, b;
	double[][] A, B;
	int numX, numY, numConstraint;
	double[] yValues; // subproblem y values
	IloRange[] subConstraint;
	HashMap<IloConstraint, IloNumExpr> rhs;
	Data data;

	public Bender(double[] c, double[] f, double[] b, double[][] A, double[][] B, Data data) throws IloException {
//		System.out.println("c= "+Arrays.toString(c));
//		System.out.println("f= "+Arrays.toString(f));
		System.out.println("b= "+Arrays.toString(b));
//		System.out.println("f= "+Arrays.toString(f));
		
		numY = c.length;
		numX = f.length;
		numConstraint = A.length;
		this.data = data;

		// this.c = Arrays.copyOf(c, c.length);
		// this.f = Arrays.copyOf(f, f.length);
		// this.b = Arrays.copyOf(b, b.length);
		this.c = c;
		this.f = f;
		this.b = b;

		this.A = A;
		this.B = B;

		rhs = new HashMap<IloConstraint, IloNumExpr>();

		// set up the master problem(which initially has no constraint)
		master = new IloCplex();
		BuildMaster();
		// master.setOut(null);
		// System.out.println("Initial master is "+master.toString());
		master.use(new BendersCallback());

		// set up the subproblem
		sub = new IloCplex();

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
			rhs.put(subConstraint[row], master.diff(b[row], master.scalProd(B[row], x)));
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
	}
	
	public void BuildMaster1() throws IloException {

		// System.out.println(data.numberOfTrucks);

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
					constraint1.addTerm(1, x[k * data.edgeSet.size() + edgeIndex]);
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
				constraint5.addTerm(data.edgeSet.get(edgeIndex).length, x[k * data.edgeSet.size() + edgeIndex]);
			}
			master.addGe(data.distanceLimit, constraint6);

		}

	}

	private class BendersCallback extends IloCplex.LazyConstraintCallback {

		protected void main() throws IloException {
			double zMaster = getValue(z);
			double[] xValue = getValues(x);
			//
			// System.out.println("In master: z= "+zMaster);
			// System.out.println("x= "+Arrays.toString(xValue));

			// set the supply constraint right-hand sides in the subproblem
			for (int row = 0; row < numConstraint; row++) {
				double temp = 0;
				for (int i = 0; i < numX; i++) {
					temp += B[row][i] * xValue[i];
				}
				subConstraint[row].setLB(b[row] - temp);
			}

			// System.out.println("Now subproblem is "+sub.toString());

			// solve the subproblem
			sub.solve();
			IloCplex.Status status = sub.getStatus();
			IloNumExpr expr = master.numExpr();

			// System.out.println("The subproblem is "+status.toString());

			if (status == IloCplex.Status.Infeasible) {
				// subproblem is infeasible -- add a feasibility cut
				// first step: get a Farkas certificate, corresponding to a dual ray
				// along which the dual is unbounded

				IloConstraint[] constraints = new IloConstraint[numConstraint];
				double[] coefficients = new double[numConstraint];
				sub.dualFarkas(constraints, coefficients);

				for (int row = 0; row < numConstraint; row++) {
					IloConstraint c = subConstraint[row];
					expr = master.sum(expr, master.prod(coefficients[row], rhs.get(c)));
				}

				// add a feasibility cut
				IloConstraint r = add(master.le(expr, 0));
				// System.out.println("\n>>> Adding feasibility cut: " + r + "\n");
				System.out.println("\n>>> Adding feasibility cut: " + "\n");
			} else if (status == IloCplex.Status.Optimal) {
				if (zMaster < sub.getObjValue() - FUZZ) {

					// the master problem surrogate variable underestimates the actual
					// flow cost -- add an optimality cut
					double[] mu = sub.getDuals(subConstraint);

					// System.out.println("mu= "+Arrays.toString(mu));

					for (int row = 0; row < numConstraint; row++) {
						expr = master.sum(expr, master.prod(mu[row], rhs.get(subConstraint[row])));
					}

					IloConstraint r = add((IloRange) master.ge(z, expr));
					// System.out.println("\n>>> Adding optimality cut: " + r + "\n");
					System.out.println("\n>>> Adding optimality cut: " + "\n");

				} else {
					System.out.println("\n>>> Accepting new incumbent with value " + getObjValue() + "\n");
					yValues = sub.getValues(y);
				}

			} else {
				// unexpected status -- report but do nothing
				System.err.println("\n!!! Unexpected subproblem solution status: " + status + "\n");
			}

		}
	}

	public final void solve() throws IloException {
//		master.setParam(IloCplex.Param.Emphasis.Memory, true);
		master.setParam(IloCplex.DoubleParam.TreLim, 2048);
		master.setParam(IloCplex.IntParam.NodeFileInd, 2);
		master.setParam(IloCplex.Param.WorkMem, 2048);
		if (master.solve()) {
			System.out.println("optimal obj= " + master.getObjValue());
			// double[] xValues = master.getValues(x);
			// System.out.println("x= " + Arrays.toString(xValues));
			// System.out.println("y= " + Arrays.toString(yValues));
		}else {
			System.out.println("The last master's status is "+master.getStatus().toString());
		}
	}

	public static void main(String[] args) throws IOException, IloException {
		Data data = new Data();
		data.readData("temp.txt");
		System.out.println("Read data done!");
//		data.readData("out_small.txt");
		data.graphTransfer();
		System.out.println("Graph transfer done!");
		data.matrixGenerator();
		System.out.println("MatrixGenerator done!");

		Bender test = new Bender(data.c, data.f, data.bb, data.A, data.B, data);
		test.solve();
	}

}
