
from fractions import Fraction
from math import lcm
import typing
from chemical_eqn import ChemEquation, Group
from matrices import Matrix
from utils import reduce_seq_to_lowest

class CouldNotSolveException(Exception):
    pass

class Solver:
    def __init__(self, equation: ChemEquation) -> None:
        self.equation = equation       
    
    def balance(self) -> ChemEquation:
        raise NotImplementedError()


def resolve_group_coefs(
    grp: Group, grp_type: typing.Literal["left", "right", "const"]
) -> dict[str, int]:
    match grp_type:
        case "const" | "right":
            return {atom: -num for atom, num in grp.atom_count_mapping.items()}
        case "left":
            return grp.atom_count_mapping



class MatrixSolver(Solver):
    def __init__(self, equation: ChemEquation) -> None:
        super().__init__(equation)
        self.const_grp, *lhs = equation.lhs
        self.rhs, self.lhs = equation.rhs, lhs

        group_count = len(self.lhs) + len(self.rhs)  # not including the constant group

        mapping = equation.atom_count_mapping()

        self.atom_dict = {atom: [0] * (group_count) for atom in mapping}
        self.const_dict = {atom: 0 for atom in mapping}
        self.progression_index = 0

    def update_dict(self, d: dict[str, int], is_const: bool = False) -> None:
        if is_const:
            self.const_dict.update(d)
        else:
            for atom, num in d.items():
                self.atom_dict[atom][self.progression_index] = num

    def balance(self) -> ChemEquation:
        self.update_dict(resolve_group_coefs(self.const_grp, "const"), is_const=True)
        for l in self.lhs:
            self.update_dict(resolve_group_coefs(l, "left"))
            self.progression_index += 1
        for r in self.rhs:
            self.update_dict(resolve_group_coefs(r, "right"))
            self.progression_index += 1

        # I am doing this to preserve order of the atoms in the dict
        order_lookup = tuple(set((*self.atom_dict.keys(), *self.const_dict.keys())))

        ab_dict: dict[tuple[float], list[int]] = {}
        for atom in order_lookup:
            const_atom = self.const_dict[atom]
            other_atom = self.atom_dict[atom]
            key = reduce_seq_to_lowest(other_atom) if const_atom == 0 else other_atom
            ab_dict[tuple(key)] = [const_atom]

        # get the coefficent groups with the least number of zeroes (most useful to us)
        unknowns_count = len(list(ab_dict.keys())[0])
        ab = sorted(
            ab_dict.items(), key=lambda eq: sum(bool(i) for i in eq[0]), reverse=True
        )[:unknowns_count] # trim off least useful coefficients and create a square matrix

        # unzipping
        a, b = Matrix([list(i[0]) for i in ab]), Matrix([i[1] for i in ab])  # type: ignore
        m = Matrix.cramers_linear_solve(a, b).transpose

        fracs = [Fraction(i).limit_denominator(15) for i in m.get_row(0)]

        l = lcm(*(i.denominator for i in fracs))
        coefs = [l, *(int(i * l) for i in fracs)]
        grps = [self.const_grp, *self.lhs, *self.rhs]

        for i in range(len(grps)):
            grps[i].coef = coefs[i]

        lhs, rhs = grps[0 : 1 + len(self.lhs)], grps[1 + len(self.lhs) :]
        final_eq = ChemEquation(lhs, rhs)
        if final_eq.is_balanced:
            return final_eq
        else:
            raise CouldNotSolveException()

