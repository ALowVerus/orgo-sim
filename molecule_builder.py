import ctypes
from config import ROOTS, RADICALS, PERIODIC_TABLE


def obj(n):
    return ctypes.cast(n, ctypes.py_object).value


def order_elements(elements):
    elts_order = []
    for elt in ["C", "H", "O"]:
        if elt in elements:
            elts_order.append(elt)
            elements.remove(elt)
    elts_order.extend(sorted(elements))
    return elts_order


def stringify_atoms_dict(formula_parts):
    elts_order = order_elements(set(formula_parts.keys()))
    formula = "".join(
        [elt + str(formula_parts[elt]) if formula_parts[elt] != 1 else elt for elt in elts_order])
    return formula


class InvalidBond(Exception):
    pass


class LockedMolecule(Exception):
    pass


class UnlockedMolecule(Exception):
    pass


class EmptyMolecule(Exception):
    pass


class Atom:

    def __init__(self, elt, idn):
        self.element = elt
        self.id = idn
        self.connections = {}

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __str__(self):
        elts_set = set([obj(n).element for n in self.connections]).difference({"H"})
        elts_order = order_elements(elts_set)
        conns_formatted = []
        for elt in elts_order:
            matches = sorted([n for n in self.connections if obj(n).element == elt], key=lambda x: obj(x).id)
            for match in matches:
                for i in range(self.connections[match]):
                    conns_formatted.append(obj(match).element + str(obj(match).id))
        conns_formatted.extend(["H" for n in self.connections if obj(n).element == "H"])
        conns = ",".join(conns_formatted)
        return "Atom({}.{}{}{})".format(self.element, self.id, ": " * (len(conns) > 0), conns)


def connect(a1, a2):
    if a1 == a2:
        raise InvalidBond
    for u, t in [(a1, a2), (a2, a1)]:
        if sum([u.connections[n] for n in u.connections]) == PERIODIC_TABLE[u.element]["valence"]:
            raise InvalidBond
    for u, t in [(a1, a2), (a2, a1)]:
        if id(t) not in u.connections:
            u.connections[id(t)] = 0
        u.connections[id(t)] += 1


class Molecule:
    def __init__(self, name=""):
        print("Initiating!")
        self.name = name
        self.atoms = []
        self.branches = []
        self.closed = False
        self._molecular_weight = 0
        self._formula = {}

    @property
    def formula(self):
        if not self.closed:
            raise UnlockedMolecule
        return self._formula

    @property
    def molecular_weight(self):
        if not self.closed:
            raise UnlockedMolecule
        return self._molecular_weight

    def create(self, element):
        if self.closed:
            raise LockedMolecule
        self.atoms.append(Atom(element, len(self.atoms) + 1))

    def brancher(self, *new_branches):
        self.log("Brancher", new_branches)
        if self.closed:
            raise LockedMolecule
        for new_branch in new_branches:
            self.create("C")
            self.branches.append([id(self.atoms[-1])])
            for i in range(new_branch - 1):
                self.create("C")
                connect(self.atoms[-1], self.atoms[-2])
                self.branches[-1].append(id(self.atoms[-1]))
        return self

    def bounder(self, *bounds):
        self.log("Bounder", bounds)
        if self.closed:
            raise LockedMolecule
        for bound in bounds:
            c1, b1, c2, b2 = bound
            connect(obj(self.branches[b1 - 1][c1 - 1]), obj(self.branches[b2 - 1][c2 - 1]))
        return self

    def mutate(self, *mutations):
        self.log("Mutating", mutations)
        if self.closed:
            raise LockedMolecule
        for mutation in mutations:
            nc, nb, elt = mutation
            target = obj(self.branches[nb - 1][nc - 1])
            if PERIODIC_TABLE[elt]["valence"] < sum([target.connections[neighbor] for neighbor in target.connections]):
                raise InvalidBond
            target.element = elt
        return self

    def add(self, *additions):
        self.log("Adding", additions)
        if self.closed:
            raise LockedMolecule
        for addition in additions:
            nc, nb, elt = addition
            target = obj(self.branches[nb - 1][nc - 1])
            if PERIODIC_TABLE[target.element]["valence"] <= sum(
                    [target.connections[neighbor] for neighbor in target.connections]):
                raise InvalidBond
            self.create(elt)
            connect(target, self.atoms[-1])
        return self

    def add_chaining(self, *args):
        self.log("Adding chains", args)
        if self.closed:
            raise LockedMolecule
        nc, nb, newbs = args[0], args[1], args[2:]
        for newb in newbs[:-1]:
            if PERIODIC_TABLE[newb]["valence"] == 1:
                raise InvalidBond
        root = obj(self.branches[nb - 1][nc - 1])
        if PERIODIC_TABLE[root.element]["valence"] <= sum(
                [root.connections[neighbor] for neighbor in root.connections]):
            raise InvalidBond
        for i, elt in enumerate(newbs):
            self.create(elt)
            connect(root, self.atoms[-1])
            root = self.atoms[-1]
        return self

    def closer(self):
        self.log("Closing")
        if self.closed:
            raise LockedMolecule
        formula_parts = {}
        for atom in self.atoms.copy():
            for i in range(
                    PERIODIC_TABLE[atom.element]["valence"] - sum([atom.connections[a] for a in atom.connections])):
                self.create("H")
                connect(atom, self.atoms[-1])
        for atom in self.atoms:
            if atom.element not in formula_parts:
                formula_parts[atom.element] = 0
            formula_parts[atom.element] += 1
        self._formula = stringify_atoms_dict(formula_parts)
        self.formula_dict = formula_parts
        self._molecular_weight = sum(formula_parts[elt] * PERIODIC_TABLE[elt]["weight"] for elt in formula_parts)
        self.closed = True
        return self

    def unlock(self):
        self.log("Opening")
        if not self.closed:
            raise UnlockedMolecule
        # Sift out empty branches
        i = len(self.branches) - 1
        while i >= 0:
            j = len(self.branches[i]) - 1
            while j >= 0:
                if obj(self.branches[i][j]).element == "H":
                    self.branches[i].pop(j)
                j -= 1
            if len(self.branches[i]) == 0:
                self.branches.pop(i)
            i -= 1
        if len(self.branches) == 0:
            raise EmptyMolecule
        # Sift out Hydrogens
        i = len(self.atoms) - 1
        while i >= 0:
            if self.atoms[i].element == "H":
                for n in self.atoms[i].connections:
                    obj(n).connections = {t: obj(n).connections[t] for t in obj(n).connections if obj(t).element != "H"}
                self.atoms.pop(i)
            i -= 1
        # Renumber
        for i, atom in enumerate(self.atoms):
            atom.id = i + 1
        # Set open
        self.closed = False
        return self

    def recursively_populate(self, node, root_branch_id=0, root_carbon_number=0):
        # Add the current node's title to the molecule
        this_branch_id = len(self.branches) + 1
        if node.title in RADICALS:
            self.brancher(RADICALS.index(node.title) + 1)
        elif node.title == "mercapto":
            self.brancher(1)
            self.mutate((1, this_branch_id, "S"))
        elif node.title in {"oxo", "al", "oyl", "oxy", "hydroxy", "ether"}:
            self.brancher(1)
            self.mutate((1, this_branch_id, "O"))
        elif node.title == "fluoro":
            self.brancher(1)
            self.mutate((1, this_branch_id, "F"))
        elif node.title == "chloro":
            self.brancher(1)
            self.mutate((1, this_branch_id, "Cl"))
        elif node.title == "bromo":
            self.brancher(1)
            self.mutate((1, this_branch_id, "Br"))
        elif node.title == "iodo":
            self.brancher(1)
            self.mutate((1, this_branch_id, "I"))
        elif node.title == "carboxy":
            self.brancher(2)
            self.mutate((2, this_branch_id, "O"))
            self.bounder((1, this_branch_id, 2, this_branch_id))
            self.add((1, this_branch_id, "O"))
        elif node.title == "formyl":
            self.brancher(2)
            self.mutate((2, this_branch_id, "O"))
            self.bounder((1, this_branch_id, 2, this_branch_id))
        elif node.title == "amido":
            self.brancher(1)
            self.mutate((1, len(self.branches), "N"))
            self.brancher(1)
            self.mutate((1, len(self.branches), "O"))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
        elif node.title == "amino" or node.title == "imino":
            self.brancher(1)
            self.mutate((1, this_branch_id, "N"))
        elif node.title == "phenyl":
            self.brancher(6)
            self.bounder((1, this_branch_id, 6, this_branch_id))
            self.bounder((1, this_branch_id, 2, this_branch_id))
            self.bounder((3, this_branch_id, 4, this_branch_id))
            self.bounder((5, this_branch_id, 6, this_branch_id))
        elif node.title == "phosphino":
            self.brancher(1)
            self.mutate((1, this_branch_id, "P"))
        elif node.title == "arsino":
            self.brancher(1)
            self.mutate((1, this_branch_id, "As"))
        elif node.title == "carbonyl":
            self.brancher(2)
            self.mutate((2, this_branch_id, "O"))
            self.bounder((1, this_branch_id, 2, this_branch_id))
        elif node.title == "oate" or node.title == "oicacid":
            self.brancher(1)
            self.mutate((1, len(self.branches), "O"))
            self.brancher(1)
            self.mutate((1, len(self.branches), "O"))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id + 1))
        # Link the new branch to the root at the correct location
        if root_branch_id and root_carbon_number:
            self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id))
            if node.title in {"oxo", "al", "imino", "oyl"}:
                self.bounder((root_carbon_number, root_branch_id, 1, this_branch_id))
        # Loop over substituents
        for locations, substituent in node.substituents:
            for location in locations:
                if substituent.title == "cyclo":
                    self.bounder((1, this_branch_id, len(self.branches[this_branch_id - 1]), this_branch_id))
                elif substituent.title in ROOTS:
                    for i in range(ROOTS.index(substituent.title)):
                        self.bounder((location, this_branch_id, location + 1, this_branch_id))
                else:
                    self.recursively_populate(substituent, this_branch_id, location)
        return self

    def log(self, *args):
        #         print("After last command:", [str(atom) for atom in self.atoms])
        #         name, content = args[0], args[1:]
        #         print(name, content)
        #         if name == "Adding":
        #             for trio in content[0]:
        #                 nc, nb, elt = trio
        #                 target = obj(self.branches[nb-1][nc-1])
        #                 print(target.element, target.id)
        pass
