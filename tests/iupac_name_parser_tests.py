from iupac_name_parser import *
from tests.test import *


Test.describe("Simple chains / impact of the number of C")

Test.it("hexane")

draw = """

CH3CH2CH2CH2CH2CH3

"""
molec = "hexane"
expected = {'C': 6, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 14}, "Wrong result for methane")
print("<COMPLETEDIN::>")

Test.it("methane")

draw = """

CH4

"""
molec = "methane"
expected = {'C': 1, 'H': 4}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 1, 'H': 4}")
test.assert_equals(ParseHer(molec).parse(), {'C': 1, 'H': 4}, "Wrong result for methane")
print("<COMPLETEDIN::>")

Test.it("ethane")

draw = """

CH3-CH3

"""
molec = "ethane"
expected = {'C': 2, 'H': 6}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 2, 'H': 6}")
test.assert_equals(ParseHer(molec).parse(), {'C': 2, 'H': 6}, "Wrong result for ethane")
print("<COMPLETEDIN::>")

Test.it("butane")

draw = """

CH3-CH2-CH2-CH3

"""
molec = "butane"
expected = {'C': 4, 'H': 10}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 4, 'H': 10}")
test.assert_equals(ParseHer(molec).parse(), {'C': 4, 'H': 10}, "Wrong result for butane")
print("<COMPLETEDIN::>")

Test.it("decane")

draw = """

CH3-CH2-CH2-CH2-CH2-CH2-CH2-CH2-CH2-CH3

"""
molec = "decane"
expected = {'C': 10, 'H': 22}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 10, 'H': 22}")
test.assert_equals(ParseHer(molec).parse(), {'C': 10, 'H': 22}, "Wrong result for decane")
print("<COMPLETEDIN::>")
print("<COMPLETEDIN::>")

Test.describe("Simple ramifications")

Test.it("REFERENCE: C8H18 (octane)")

draw = """

CH3-CH2-CH2-CH2-CH2-CH2-CH2-CH3

"""
molec = "octane"
expected = {'C': 8, 'H': 18}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 18}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 18}, "Wrong result for octane")
print("<COMPLETEDIN::>")

Test.it("One ramification")

draw = """

        CH2-CH3
        |
CH3-CH2-CH-CH2-CH2-CH3

"""
molec = "3-ethylhexane"
expected = {'C': 8, 'H': 18}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 18}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 18}, "Wrong result for 3-ethylhexane")
print("<COMPLETEDIN::>")

Test.it("Two ramifications")

draw = """

   CH3  CH2-CH3
    |   |
CH3-CH-CH-CH2-CH3

"""
molec = "3-ethyl-2-methylpentane"
expected = {'C': 8, 'H': 18}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 18}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 18}, "Wrong result for 3-ethyl-2-methylpentane")
print("<COMPLETEDIN::>")

Test.it("Two ramifications on the same C")

draw = """

        CH2-CH3
        |
CH3-CH2-C-CH2-CH3
        |
        CH3

"""
molec = "3-ethyl-3-methylpentane"
expected = {'C': 8, 'H': 18}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 18}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 18}, "Wrong result for 3-ethyl-3-methylpentane")
print("<COMPLETEDIN::>")

Test.it("Handle multipliers")

draw = """

        CH3
        |
CH3-CH2-C-CH2-CH2-CH3
        |
        CH3

"""
molec = "3,3-dimethylhexane"
expected = {'C': 8, 'H': 18}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 18}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 18}, "Wrong result for 3,3-dimethylhexane")
print("<COMPLETEDIN::>")

Test.it("Handle multipliers")

draw = """

 CH3   CH3
   \\   /
CH3-C-C-CH3
   /   \\
 CH3   CH3

"""
molec = "2,2,3,3-tetramethylbutane"
expected = {'C': 8, 'H': 18}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 18}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 18}, "Wrong result for 2,2,3,3-tetramethylbutane")
print("<COMPLETEDIN::>")
print("<COMPLETEDIN::>")

Test.describe("Effect of cycles and multiple bonds")

Test.it("REFERENCE: C8H16 (cyclooctane)")

draw = """

 CH2-CH2-CH2-CH2
 |           |
 CH2-CH2-CH2-CH2

"""
molec = "cyclooctane"
expected = {'C': 8, 'H': 16}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 16}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 16}, "Wrong result for cyclooctane")
print("<COMPLETEDIN::>")

Test.it("One cycle of size 6 and one ramification")

draw = """

 CH2-CH2-CH-CH2-CH3
 |       |
 CH2-CH2-CH2

"""
molec = "1-ethylcyclohexane"
expected = {'C': 8, 'H': 16}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 16}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 16}, "Wrong result for 1-ethylcyclohexane")
print("<COMPLETEDIN::>")

Test.it("One cycle of size 4 and several ramifications")

draw = """

 CH2-CH-CH2-CH3
 |   |
 CH2-C-CH3
     |
     CH3

"""
molec = "1-ethyl-2,2-dimethylcyclobutane"
expected = {'C': 8, 'H': 16}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 16}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 16}, "Wrong result for 1-ethyl-2,2-dimethylcyclobutane")
print("<COMPLETEDIN::>")

Test.it("One double bond: at an extremity")

draw = """

CH2=CH-CH2-CH2-CH2-CH2-CH2-CH3

"""
molec = "oct-1-ene"
expected = {'C': 8, 'H': 16}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 16}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 16}, "Wrong result for oct-1-ene")
print("<COMPLETEDIN::>")

Test.it("One double bond: anywhere in the chain")

draw = """

CH3-CH2-CH=CH-CH2-CH2-CH2-CH3

"""
molec = "oct-3-ene"
expected = {'C': 8, 'H': 16}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 16}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 16}, "Wrong result for oct-3-ene")
print("<COMPLETEDIN::>")

Test.it("One double bond: elision of the position '-1-'")

draw = """

CH2=CH-CH2-CH2-CH2-CH2-CH2-CH3

"""
molec = "octene"
expected = {'C': 8, 'H': 16}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 16}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 16}, "Wrong result for octene")
print("<COMPLETEDIN::>")
print("<COMPLETEDIN::>")

Test.describe("Effect of mutliple bonds and cycles, part 2")

Test.it("Double bonds")

draw = """

CH3-CH=CH-CH2-CH=CH-CH2-CH3

"""
molec = "oct-2,5-diene"
expected = {'C': 8, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 14}, "Wrong result for oct-2,5-diene")
print("<COMPLETEDIN::>")

Test.it("Triple bond: at an extremity")

draw = """

CH{=}C-CH2-CH2-CH2-CH2-CH2-CH3      "{=}" used as triple bond (should be 3 lines)

"""
molec = "oct-1-yne"
expected = {'C': 8, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 14}, "Wrong result for oct-1-yne")
print("<COMPLETEDIN::>")

Test.it("Triple bond: in the chain")

draw = """

CH3-C{=}C-CH2-CH2-CH2-CH2-CH3

"""
molec = "oct-2-yne"
expected = {'C': 8, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 14}, "Wrong result for oct-2-yne")
print("<COMPLETEDIN::>")

Test.it("Triple bond: elision of the position")

draw = """

CH{=}C-CH2-CH2-CH2-CH2-CH2-CH3

"""
molec = "octyne"
expected = {'C': 8, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 14}, "Wrong result for octyne")
print("<COMPLETEDIN::>")

Test.it("Mix of cycles and multiple bonds")

draw = """

 CH2-CH2-CH-CH2-CH3
 |       |
 CH=CH-CH2

"""
molec = "3-ethylcyclohexene"
expected = {'C': 8, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 8, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 8, 'H': 14}, "Wrong result for 3-ethylcyclohexene")
print("<COMPLETEDIN::>")
print("<COMPLETEDIN::>")

Test.describe("Simple functions: oxygen")

Test.it("REFERENCE: C5H12 (pentane)")

draw = """

CH3-CH2-CH2-CH2-CH3

"""
molec = "pentane"
expected = {'C': 5, 'H': 12}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 12}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 12}, "Wrong result for pentane")
print("<COMPLETEDIN::>")

Test.it("pentanol")

draw = """

CH3-CH2-CH2-CH2-CH2-OH

"""
molec = "pentanol"
expected = {'C': 5, 'H': 12, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 12, 'O': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 12, 'O': 1}, "Wrong result for pentanol")
print("<COMPLETEDIN::>")

Test.it("pentan-2-ol")

draw = """

    OH
    |
CH3-CH-CH2-CH2-CH3

"""
molec = "pentan-2-ol"
expected = {'C': 5, 'H': 12, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 12, 'O': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 12, 'O': 1}, "Wrong result for pentan-2-ol")
print("<COMPLETEDIN::>")

Test.it("pentan-2,4-diol")

draw = """

    OH     OH
    |      |
CH3-CH-CH2-CH-CH3

"""
molec = "pentan-2,4-diol"
expected = {'C': 5, 'H': 12, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 12, 'O': 2}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 12, 'O': 2}, "Wrong result for pentan-2,4-diol")
print("<COMPLETEDIN::>")

Test.it("pentanal")

draw = """

CH3-CH2-CH2-CH2-CH=O

"""
molec = "pentanal"
expected = {'C': 5, 'H': 10, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 10, 'O': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 10, 'O': 1}, "Wrong result for pentanal")
print("<COMPLETEDIN::>")

Test.it("pentan-2-one")

draw = """

    O
    ||
CH3-C-CH2-CH2-CH3

"""
molec = "pentan-2-one"
expected = {'C': 5, 'H': 10, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 10, 'O': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 10, 'O': 1}, "Wrong result for pentan-2-one")
print("<COMPLETEDIN::>")

Test.it("pentandial")

draw = """

O=CH-CH2-CH2-CH2-CH=O

"""
molec = "pentandial"
expected = {'C': 5, 'H': 8, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 8, 'O': 2}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 8, 'O': 2}, "Wrong result for pentandial")
print("<COMPLETEDIN::>")

Test.it("pentan-2,4-dione")

draw = """

    O     O
    ||    ||
CH3-C-CH2-C-CH3

"""
molec = "pentan-2,4-dione"
expected = {'C': 5, 'H': 8, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 8, 'O': 2}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 8, 'O': 2}, "Wrong result for pentan-2,4-dione")
print("<COMPLETEDIN::>")
print("<COMPLETEDIN::>")

Test.describe("Simple functions: halogens")

Test.it("REFERENCE: C5H12 (pentane)")

draw = """

CH3-CH2-CH2-CH2-CH3

"""
molec = "pentane"
expected = {'C': 5, 'H': 12}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'H': 12}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'H': 12}, "Wrong result for pentane")
print("<COMPLETEDIN::>")

Test.it("1-fluoropentane")

draw = """

CH3-CH2-CH2-CH2-CH2-F

"""
molec = "1-fluoropentane"
expected = {'C': 5, 'F': 1, 'H': 11}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'F': 1, 'H': 11}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'F': 1, 'H': 11}, "Wrong result for 1-fluoropentane")
print("<COMPLETEDIN::>")

Test.it("2-chloropentane")

draw = """

    Cl
    |
CH3-CH-CH2-CH2-CH3

"""
molec = "2-chloropentane"
expected = {'C': 5, 'Cl': 1, 'H': 11}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 5, 'Cl': 1, 'H': 11}")
test.assert_equals(ParseHer(molec).parse(), {'C': 5, 'Cl': 1, 'H': 11}, "Wrong result for 2-chloropentane")
print("<COMPLETEDIN::>")

Test.it("1-bromo-4-chloropentane")

draw = """

    Cl
    |
CH3-CH-CH2-CH2-CH2-Br

"""
molec = "1-bromo-4-chloropentane"
expected = {'Br': 1, 'C': 5, 'Cl': 1, 'H': 10}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'Br': 1, 'C': 5, 'Cl': 1, 'H': 10}")
test.assert_equals(ParseHer(molec).parse(), {'Br': 1, 'C': 5, 'Cl': 1, 'H': 10},
                   "Wrong result for 1-bromo-4-chloropentane")
print("<COMPLETEDIN::>")
print("<COMPLETEDIN::>")

Test.describe("Simple functions: nitrogen")

Test.it("REFERENCE: C6H14 (hexane)")

draw = """

CH3-CH2-CH2-CH2-CH2-CH3

"""
molec = "hexane"
expected = {'C': 6, 'H': 14}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 14}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 14}, "Wrong result for hexane")
print("<COMPLETEDIN::>")

Test.it("hexylamine")

draw = """

CH3-CH2-CH2-CH2-CH2-CH2-NH2

"""
molec = "hexylamine"
expected = {'C': 6, 'H': 15, 'N': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 15, 'N': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 15, 'N': 1}, "Wrong result for hexylamine")
print("<COMPLETEDIN::>")

Test.it("butylethylamine")

draw = """

CH3-CH2-CH2-CH2-NH-CH2-CH3

"""
molec = "butylethylamine"
expected = {'C': 6, 'H': 15, 'N': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 15, 'N': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 15, 'N': 1}, "Wrong result for butylethylamine")
print("<COMPLETEDIN::>")

Test.it("ethylmethylpropylamine")

draw = """

            CH3
            |
CH3-CH2-CH2-N-CH2-CH3

"""
molec = "ethylmethylpropylamine"
expected = {'C': 6, 'H': 15, 'N': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 15, 'N': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 15, 'N': 1}, "Wrong result for ethylmethylpropylamine")
print("<COMPLETEDIN::>")

Test.it("triethylamine")

draw = """

N(CH2-CH3)3

"""
molec = "triethylamine"
expected = {'C': 6, 'H': 15, 'N': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 15, 'N': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 15, 'N': 1}, "Wrong result for triethylamine")
print("<COMPLETEDIN::>")

Test.it("Alternative nomenclature: hexan-1,6-diamine")

draw = """

NH2-CH2-CH2-CH2-CH2-CH2-CH2-NH2

"""
molec = "hexan-1,6-diamine"
expected = {'C': 6, 'H': 16, 'N': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 16, 'N': 2}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 16, 'N': 2}, "Wrong result for hexan-1,6-diamine")
print("<COMPLETEDIN::>")

Test.it("WARNING: amiDe, not amiNe, here!")

draw = """

                    O
                    ||
CH3-CH2-CH2-CH2-CH2-C-NH2

"""
molec = "hexanamide"
expected = {'C': 6, 'H': 13, 'N': 1, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 6, 'H': 13, 'N': 1, 'O': 1}")
test.assert_equals(ParseHer(molec).parse(), {'C': 6, 'H': 13, 'N': 1, 'O': 1}, "Wrong result for hexanamide")
print("<COMPLETEDIN::>")

Test.it("REFERENCE: C14H30O2 (tetradecandiol)")

draw = """

CH2-CH2-CH2-CH2-CH2-CH2-CH2-OH
|
CH2-CH2-CH2-CH2-CH2-CH2-CH2-OH

"""
molec = "tetradecandiol"
expected = {'C': 14, 'H': 30, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + "{'C': 14, 'H': 30, 'O': 2}")
test.assert_equals(ParseHer(molec).parse(), {'C': 14, 'H': 30, 'O': 2}, "Wrong result for hexanamide")
print("<COMPLETEDIN::>")

Test.it("4-[1-methyl]ethyl-6-hydroxyhex-2,4-dienoic acid")

draw = """

   O   CH3-CH-CH3
  ||       |
HO-C-CH=CH-C=CH-CH2-OH

"""
molec = "4-[1-methyl]ethyl-6-hydroxyhex-2,4-dienoic acid"
expected = {'C': 9, 'H': 14, 'O': 3}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")

Test.it("propanoic acid")

draw = """

        O
        ||
CH3-CH2-C-OH

"""
molec = "propanoic acid"
expected = {'C': 3, 'H': 6, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")

Test.it("hex-2,4-dienoic acid")

draw = """

        O
        ||
CH3-CH2-C-OH

"""
molec = "hex-2,4-dienoic acid"
expected = {'C': 6, 'H': 8, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")

Test.it("propdial")

draw = """

CHO-CH2-CH3

"""
molec = "propal"
expected = {'C': 3, 'H': 6, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")

Test.it("propdial")

draw = """

CHO-CH2-CHO

"""
molec = "propdial"
expected = {'C': 3, 'H': 4, 'O': 2}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")

Test.it("bromo-5,10,13-tridecyltetradec-1,13-diene")

draw = """

Br            C10H21             C10H21     C10H21
|             |                  |          |
CH=CH-CH2-CH2-CH-CH2-CH2-CH2-CH2-CH-CH2-CH2-C=CH2

"""
molec = "bromo-5,10,13-tridecyltetradec-1,13-diene"
expected = {'C': 44, 'H': 85, 'Br': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")

Test.it("bromo-5-tridecyltetradec-1,13-diene")

draw = """

Br            C13H27
|             |
CH=CH-CH2-CH2-CH-CH2-CH2-CH2-CH2-CH2-CH2-CH2-CH=CH2

"""
molec = "bromo-5-tridecyltetradec-1,13-diene"
expected = {'H': 51, 'C': 27, 'Br': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")


Test.it("BLOOOOP")

draw = """
BLOOP
"""
molec = "2-hexadec-15-en-1,3,5,12,13-pentaynyl-8-hexadecyl-2,3,4,4,9-pentanon-5-enyl-8-prop-2-ynylcyclodec-1,5,6,7-tetrayne"
expected = {'C': 90, 'H': 128}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


Test.it("Methoxies")

draw = """
BLOOP
"""
molec = "7-hydroxy-2-methoxycyclononanone"
expected = {'C': 10, 'H': 18, "O": 3}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


Test.it("-oxycarbonyl-")

draw = """
BLOOP
"""
molec = "3-[4-carboxy-2-methoxycarbonyl]butyldodecandioic acid"
expected = {'C': 19, 'H': 32, "O": 8}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


Test.it("triphenylarsine")

draw = """
BLOOP
"""
molec = "triphenylarsine"
expected = {'C': 18, 'H': 15, "As": 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


name = "3,4,5,8,10,14-hexa[2,4,7-tri[11-hexadecyl]tridec-2,6,7,11,12-pentaenyl]cyclodecylcycloheptadec-2,9-diene"
Test.it(name)

draw = """
C
"""
molec = name
expected = {'C': 599, 'H': 1002}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result.")
print("<COMPLETEDIN::>")


Test.it("weird elisions")

draw = """
BLOOP
"""
molec = "2-[4,6,9,11-tetra[[10-phenyl-11-[[6,16,18-triarsino]octadecyl[2-iodo]ethyl]phosphino]tetradecyl]amino-10-chloro]tetradecoxycarbonylbutandioic acid"
expected = {'C': 179, 'H': 337, 'As': 12, 'I': 4, 'P': 4, 'N': 4, 'Cl': 1, 'O': 6}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


Test.it("-oate")

draw = """
BLOOP
"""
molec = "cyclobutyl 6-oxoheptanoate"
expected = {'C': 11, 'H': 18, "O": 3}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


Test.it("-dioate")

draw = """
BLOOP
"""
molec = "butyl ethan-1,2-dioate"
expected = {'C': 10, 'H': 18, 'O': 4}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")


Test.it("ether")

draw = """
BLOOP
"""
molec = "[6-[2-hexyl-1,3,7-triiodo]heptyl]decyl[6,9-diphosphino]undecylether"
expected = {'C': 34, 'H': 69, 'I': 3, 'P': 2, 'O': 1}

print(molec + ":\n" + draw.strip("\n") + "\n" + str(expected))
test.assert_equals(ParseHer(molec).parse(), expected, "Wrong result for BLOOP")
print("<COMPLETEDIN::>")

