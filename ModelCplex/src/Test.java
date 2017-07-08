import java.util.Scanner;
import java.util.regex.MatchResult;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.cplex.IloCplex;

public class Test {

	public static void main(String[] args) throws IloException {
		// String input = "1 fish 2 fish red fish blue fish";
		// Scanner s = new Scanner(input);
		// String out=s.findInLine("(\\d+) fish (\\d+) fish (\\w+) fish
		// (\\w+)");
		// MatchResult result = s.match();
		// for (int i=1; i<=result.groupCount(); i++)
		// System.out.println(result.group(i));
		// s.close();
		// System.out.println(out.toString());
		// String s=" 8";
		// System.out.println(Double.valueOf(s));

		IloCplex cplex = new IloCplex();
		IloLinearNumExpr constraint = cplex.linearNumExpr();
		IloIntVar[] x = new IloIntVar[3];
		x[0]=cplex.intVar(0, 1,"x1");
		
		constraint=cplex.linearNumExpr();
		constraint.addTerm(1, x[0]);
		System.out.println(constraint.toString());
	}
}
