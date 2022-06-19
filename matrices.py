import typing
from dataclasses import dataclass

num = int | float


@dataclass
class Point:
    x: num
    y: num


class Matrix:
    @staticmethod
    def _to_matrix(m: typing.Union["Matrix", list[list[num]]]) -> "Matrix":
        return m if isinstance(m, Matrix) else Matrix(m)

    def __init__(
        self, data: list[list[num]] | None = None, n: int | None = None
    ) -> None:
        if data:
            if n is not None and any(len(i) != n for i in data):
                raise AttributeError("number of columns is not equal to n")
            self._n = len(data[0])
            # we don't make a self._m because that can be calculated
            self._data = data
        else:
            if n is None:
                raise AttributeError("n is required when data is not provided")
            else:
                self._n = n
                self._data = []

    def resolve_coord(self, val: int | tuple[num, num] | None) -> Point:
        match val:
            case int(x):
                return Point(0, x)
            case (x, y):
                return Point(x, y)
            case None:
                return Point(self.n - 1, self.m - 1)

    def all_except(self, p: Point) -> "Matrix":
        return Matrix(
            [
                [x for i, x in enumerate(y) if i != p.x]
                for j, y in enumerate(self._data)
                if j != p.y
            ]
        )

    def minor(self, p: Point) -> num:
        return self.all_except(p).determinant

    def cofactor(self, p: Point) -> num:
        x, y = int(p.x or 0), int(p.y or 0)  # strict type checking moment
        if x is not None and y is not None:
            return ((-1) ** (x + y)) * self.minor(Point(x, y))
        else:
            raise ValueError("no None coords")

    @property
    def m(self) -> int:
        return len(self._data)

    @property
    def n(self) -> int:
        return len(self._data[0]) if self._data else self._n

    @property
    def determinant(self) -> num:
        if self.m != self.n:
            raise ValueError(
                "Determinants can only be calculated with square matrices (m=n)"
            )
        if self.m == 1:
            return self._data[0][0]
        return sum(
            ((-1) ** x) * self.get((x, 0)) * self.all_except(Point(x, 0)).determinant
            for x in range(self.n)
        )

    @property
    def minor_matrix(self) -> "Matrix":
        return Matrix(
            [[self.minor(Point(i, j)) for i in range(self.m)] for j in range(self.m)]
        )

    @property
    def cofactor_matrix(self) -> "Matrix":
        return Matrix(
            [[self.cofactor(Point(i, j)) for i in range(self.m)] for j in range(self.m)]
        )

    @staticmethod
    def identity(n: int) -> "Matrix":
        return Matrix([[1 if i == j else 0 for j in range(n)] for i in range(n)])

    @property
    def transpose(self) -> "Matrix":
        return Matrix(
            [[self._data[i][j] for i in range(self.m)] for j in range(self.n)]
        )

    @property
    def adjoint(self) -> "Matrix":
        return self.cofactor_matrix.transpose

    @property
    def inverse(self) -> "Matrix":
        return self.adjoint / self.determinant

    @staticmethod
    def cramers_linear_solve(a: "Matrix", b: "Matrix") -> "Matrix":
        d = a.determinant
        ds: list[Matrix] = []
        for i in range(a.n):
            m = a.copy()
            m[(i, 0) : (i + 1, a.m)] = b
            ds.append(m)
        return Matrix([[i.determinant / d for i in ds]]).transpose

    def get(self, coord: tuple[int, int]) -> num:
        return self._data[coord[1]][coord[0]]
    
    def get_row(self, m: int) -> list[num]:
        return self._data[m]

    def copy(self) -> "Matrix":
        return Matrix([i.copy() for i in self._data])

    def __mul__(self, other: num) -> "Matrix":
        return Matrix([[other * i for i in j] for j in self._data])

    def __truediv__(self, other: num) -> "Matrix":
        return Matrix([[i / other for i in j] for j in self._data])

    def __str__(self) -> str:
        return "\n".join(" ".join(str(j) for j in i) for i in self._data)

    def __getitem__(
        self, val: slice 
    ) -> "Matrix":
        start, end = self.resolve_coord(val.start), self.resolve_coord(val.stop)
        return Matrix([i[start.x : end.x] for i in self._data[start.y : end.y]])

    def __setitem__(
        self, key: slice, value: typing.Union["Matrix", list[list[num]]]
    ) -> None:
        matrix = Matrix._to_matrix(value)
        start, stop = self.resolve_coord(key.start), self.resolve_coord(key.stop)
        width, height = stop.x - start.x, stop.y - start.y
        if matrix.n == width and matrix.m == height:
            for j in range(matrix.m):
                for i in range(matrix.n):
                    self._data[int(start.y) + j][int(start.x) + i] = matrix.get((i, j))
        else:
            raise ValueError("Cannot set value of unequal dimensions")


if __name__ == "__main__":
    x = Matrix(
        [
            [1, 2, 3],
            [1, 5, 5],
        ]
    )

    x[(1, 0):] = Matrix(
        [
            [4, 4],
            [5, 5],
        ]
    )
    print(x)
