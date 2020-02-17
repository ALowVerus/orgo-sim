
# FOR REFERENCE
RADICALS = ["meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec", "undec", "dodec",
            "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec", "nonadec"]
MULTIPLIERS = ["", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca", "undeca", "dodeca",
               "trideca", "tetradeca", "pentadeca", "hexadeca", "heptadeca", "octadeca", "nonadeca"]
ROOTS = [
    "an", "en", "yn"
]
REDUNDANT_SUFFIXES = {
    "ol": "hydroxy",
    "one": "oxo",
    "amide": "amido",
    "amine": "amino",
    "imine": "imino",
    "benzene": "phenyl",
    "thiol": "mercapto",
    "phosphine": "phosphino",
    "arsine": "arsino",
    "carboxylicacid": "carboxy"
}
SUFFIXES = {
    "al": {"O": 1, "H": -2},
    "oicacid": {"O": 2, "H": -2},
    "ether": {"O": 1},
}
PREFIXES = {
    "cyclo": {"H": -2},
    "hydroxy": {"O": 1},
    "oxo": {"O": 1, "H": -2},
    "carboxy": {"O": 2, "C": 1},
    "formyl": {"C": 1, "O": 1},
    "amido": {"N": 1, "O": 1, "H": -1},
    "amino": {"N": 1, "H": 1},
    "imino": {"N": 1, "H": -1},
    "phenyl": {"C": 6, "H": 4},
    "mercapto": {"S": 1},
    "phosphino": {"P": 1, "H": 1},
    "arsino": {"As": 1, "H": 1},
    "fluoro": {"F": 1, "H": -1},
    "chloro": {"Cl": 1, "H": -1},
    "bromo": {"Br": 1, "H": -1},
    "iodo": {"I": 1, "H": -1}
}
MEDIXES = {
    "oxy": {"O": 1},
    "yl": {},
    "oyl": {"O": 1, "H": -2},
    "carbonyl": {"C": 1, "O": 1},
    "oate": {"O": 2, "H": -2},
}
PERIODIC_TABLE = {
    "H": {"valence": 1, "weight": 1},
    "B": {"valence": 3, "weight": 10.8},
    "C": {"valence": 4, "weight": 12},
    "N": {"valence": 3, "weight": 14},
    "As": {"valence": 3, "weight": 74.9},
    "O": {"valence": 2, "weight": 16},
    "F": {"valence": 1, "weight": 19},
    "Mg": {"valence": 2, "weight": 24.3},
    "P": {"valence": 3, "weight": 31},
    "S": {"valence": 2, "weight": 32.1},
    "Cl": {"valence": 1, "weight": 35.5},
    "Br": {"valence": 1, "weight": 80},
    "I": {"valence": 1, "weight": 126.9}
}