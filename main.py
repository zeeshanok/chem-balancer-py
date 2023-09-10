from chemical_eqn import ChemEquation
from solvers import CouldNotSolveException, MatrixSolver

CURSOR_UP = "\x1B[1A"
CLEAR_LINE = "\x1B[2K"

if __name__ == "__main__":
    while True:
        try:
            rxn = ChemEquation.parse(input("> "))
            print(f"{CURSOR_UP}{CLEAR_LINE}In: {rxn}")
            solver = MatrixSolver(rxn)
            print(f"Out: {solver.balance()}")
            print()
        except CouldNotSolveException as e:
            print(f"We couldn't balance your equation\n{e}")   
        except Exception as e:
            print(e)
