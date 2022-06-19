from fractions import Fraction
from math import lcm
import typing
from chemical_eqn import ChemEquation, Group
from matrices import Matrix, num

def resolve_group_coefs(
    grp: Group, grp_type: typing.Literal["left", "right", "const"]
) -> dict[str, int]:
    match grp_type:
        case "const" | "right":
            return {atom: -num for atom, num in grp.atom_count_mapping.items()}
        case "left":
            return grp.atom_count_mapping


class Solver:
    def __init__(self, eqn: ChemEquation) -> None:
        self.const_grp, *lhs = eqn.lhs
        self.rhs, self.lhs = eqn.rhs, lhs
        group_count = len(self.lhs) + len(self.rhs) + 1
        if eqn.atom_count != group_count - 1:
            raise ValueError("This equation cannot be balanced (yet)")

        mapping = eqn.atom_count_mapping()
        self.atom_dict = {atom: [0] * (group_count - 1) for atom in mapping}
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

        a_l: list[list[num]] = []
        b_l: list[list[num]] = []
        for atom in order_lookup:
            a_l.append([int(i) for i in self.atom_dict[atom]])
            b_l.append([self.const_dict[atom]])
        
        a, b = Matrix(a_l), Matrix(b_l)
        m = Matrix.cramers_linear_solve(a, b).transpose

        fracs = [Fraction(i).limit_denominator(15) for i in  m.get_row(0)]
        l = lcm(*(i.denominator for i in fracs))
        coefs = [l, *(int(i * l) for i in fracs)]
        grps = [self.const_grp, *self.lhs, *self.rhs]
        for i in range(len(grps)):
            grps[i].coef = coefs[i]
       
        lhs, rhs = grps[0:1+len(self.lhs)], grps[1+len(self.lhs):]
        return ChemEquation(lhs, rhs)



if __name__ == "__main__":
    while True:
        try:
            rxn = ChemEquation.parse(input("> "))
            print(rxn)
            solver = Solver(rxn)
            print(solver.balance())
        except Exception as e:
            print(e)
