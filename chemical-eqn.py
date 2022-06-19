from collections import Counter
from dataclasses import dataclass
from functools import reduce
from typing import Callable, Union

ZERO_SUB = ord("₀")


def get_subscript(n: int) -> str:
    return "".join(chr(ZERO_SUB + int(i)) for i in str(n))


def get_closing_bracket(b: str) -> str | None:
    return {
        "(": ")",
        "[": "]",
        "{": "}",
    }.get(b)


def coef_convert(v: int, n: int = 1) -> str:
    return str(v) if v > n else ""


def subscript_coef_convert(v: int, n: int = 1) -> str:
    return get_subscript(v) if v > n else ""


@dataclass
class ElementGroup:
    coef: int
    atom: str

    def __repr__(self) -> str:
        return f"{self.atom}{subscript_coef_convert(self.coef)}"

    def __mul__(self, other: int) -> "ElementGroup":
        return ElementGroup(self.coef * other, self.atom)


class Group:
    def __init__(
        self,
        coef: int,
        children: list[Union[ElementGroup, "Group"]],
        is_sub: bool = False,
    ) -> None:
        self.coef = coef
        self.is_sub = is_sub
        if children:
            self.children = children
        else:
            raise ValueError("Children cannot be empty")

    @property
    def atom_count_mapping(self) -> dict[str, int]:
        """Returns a dict mapping every atom type to its number"""
        counter: Counter[str] = Counter()
        for child in self.children:
            match child:
                case ElementGroup() as g:
                    counter[g.atom] += g.coef
                case Group() as g:
                    counter.update(g.atom_count_mapping)
        return {atom: self.coef * coef for atom, coef in counter.items()}

    @property
    def atom_count(self) -> int:
        """Returns number of atom types in this group"""
        return len(self.atom_count_mapping)

    @staticmethod
    def parse(s: str) -> "Group":
        lexer = GroupLexer(s)
        return lexer.get_group()

    def __str__(self) -> str:
        s = "".join(str(i) for i in self.children)
        return (
            f"({s}){subscript_coef_convert(self.coef)}"
            if self.is_sub
            else f"{coef_convert(self.coef)}{s}"
        )

    def __repr__(self) -> str:
        return self.__str__()


@dataclass
class ChemEquation:
    lhs: list[Group]
    rhs: list[Group]

    def atom_count_mapping(self, is_left: bool) -> dict[str, int]:
        return reduce(
            lambda x, y: x + Counter(y),
            (i.atom_count_mapping for i in (self.lhs if is_left else self.rhs)),
            Counter[str](),
        )
    
    @property
    def is_balanced(self) -> bool:
        '''Returns `True` if the chemical equation is balanced otherwise `False`'''
        return self.atom_count_mapping(True) == self.atom_count_mapping(False)

    @staticmethod
    def parse(s: str) -> "ChemEquation":
        l, r, *_ = s.split("->")
        return ChemEquation(
            [Group.parse(i.strip()) for i in l.split("+")],
            [Group.parse(i.strip()) for i in r.split("+")],
        )

    def __repr__(self) -> str:
        l = " + ".join(str(i) for i in self.lhs)
        r = " + ".join(str(i) for i in self.rhs)
        return f"{l} -> {r}"


class GroupLexer:
    def __init__(self, text: str) -> None:
        self.text = text
        self._index = 0

    @property
    def at_end(self) -> bool:
        return self._index > len(self.text) - 1

    def consume(self) -> str:
        self._index += 1
        return self.text[self._index - 1]

    def peek(self, n: int = 1) -> str:
        return self.text[self._index : self._index + n]

    def consume_numeric(self) -> int:
        coef_s = ""
        while self.peek().isnumeric():
            coef_s += self.consume()
        return int(coef_s)

    def consume_while(self, fn: Callable[[int, str], bool]) -> str:
        c = ""
        while not self.at_end and fn(self._index, self.peek()):
            c += self.consume()
        return c

    def consume_coef(self) -> int:
        """Returns the next number in `text` if it exists otherwise returns 1"""
        return self.consume_numeric() if self.peek().isnumeric() else 1

    def get_group(self) -> Group:
        coef = self.consume_coef()

        elems: list[ElementGroup | Group] = []

        while not self.at_end:
            c = self.consume()
            if c.isupper():
                while not self.at_end and self.peek().islower():
                    c += self.consume()
                mul = self.consume_coef()
                try:
                    elems.append(ElementGroup(mul, c))
                except:
                    continue
            if (closing_bracket := get_closing_bracket(c)) is not None:
                sub_group_text = self.consume_while(lambda _, ch: ch != closing_bracket)
                self.consume()  # consume closing bracket
                lexer = GroupLexer(sub_group_text)
                sub_group = lexer.get_group()
                sub_group.coef = self.consume_coef()
                sub_group.is_sub = True
                elems.append(sub_group)

        return Group(coef, elems)


while True:
    eqn = ChemEquation.parse(input("> "))
    print()
    print(
        f"LHS\n{eqn.atom_count_mapping(True)}\n\nRHS\n{eqn.atom_count_mapping(False)}"
    )
    print(f"This equation is{' not' if not eqn.is_balanced else ''} balanced")
