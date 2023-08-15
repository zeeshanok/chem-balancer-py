from chemical_eqn import ChemEquation
from solvers import CouldNotSolveException, MatrixSolver

if __name__ == "__main__":
    while True:
        try:
            rxn = ChemEquation.parse(input("> "))
            print(rxn)
            solver = MatrixSolver(rxn)
            print(solver.balance())
        except CouldNotSolveException as e:
            print(f"We couldn't balance your equation\n{e}")   
        except Exception as e:
            print(e)
