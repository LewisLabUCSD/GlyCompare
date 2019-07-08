import unittest

import glypy
from glypy.composition import composition_transform
from glypy.structure import glycan_composition
from glypy import monosaccharides, Substituent, glycans

water_mass = glypy.Composition("H2O").mass

GlycanComposition = glycan_composition.GlycanComposition
MonosaccharideResidue = glycan_composition.MonosaccharideResidue

FrozenGlycanComposition = glycan_composition.FrozenGlycanComposition
FrozenMonosaccharideResidue = glycan_composition.FrozenMonosaccharideResidue


class IUPACLiteTests(unittest.TestCase):
    def test_parse_identity(self):
        items = ["Glc2NAc", "Neu5Ac", "Fuc", "6-dHex"]
        for item in items:
            term = glycan_composition.from_iupac_lite(item)
            deparse = glycan_composition.to_iupac_lite(term)
            self.assertEqual(item, deparse)


class MonosaccharideResidueTests(unittest.TestCase):
    def test_from_monosaccharide(self):
        for base in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[base]
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)
            self.assertAlmostEqual(base.mass() - water_mass, residue.mass(), 3)

    def test_clone_coherence(self):
        for base in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[base]
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)
            dup = residue.clone()
            self.assertAlmostEqual(dup.mass(), residue.mass(), 3)

    def test_derivatize(self):
        reference = {"GlcNAc": 245.1263, "Fuc": 174.0892, "NeuAc": 361.1737}
        for name in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[name]
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)
            composition_transform.derivatize(base, "methyl")
            composition_transform.derivatize(residue, "methyl")
            self.assertAlmostEqual(
                base.mass() - (water_mass + 2 * (Substituent("methyl").mass() - glypy.Composition("H").mass * 2)),
                residue.mass(), 3)
            self.assertAlmostEqual(reference[name], residue.mass(), 3)

            base = monosaccharides[name]
            composition_transform.derivatize(base, "methyl")
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)

            self.assertAlmostEqual(
                base.mass() - (water_mass + 2 * (Substituent("methyl").mass() - glypy.Composition("H").mass * 2)),
                residue.mass(), 3)

            self.assertAlmostEqual(reference[name], residue.mass(), 3)

    def test_acidic_special_case(self):
        acidic_hexose1 = glycan_composition.MonosaccharideResidue.from_iupac_lite("HexA")
        acidic_hexose2 = glycan_composition.MonosaccharideResidue.from_iupac_lite("aHex")
        hexose = glycan_composition.MonosaccharideResidue.from_iupac_lite("Hex")
        self.assertAlmostEqual(acidic_hexose1.mass() - hexose.mass(), 13.9792, 3)
        self.assertAlmostEqual(acidic_hexose2.mass() - hexose.mass(), 13.9792, 3)


class FrozenMonosaccharideResidueTests(unittest.TestCase):
    def test_from_monosaccharide(self):
        for base in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[base]
            residue = glycan_composition.FrozenMonosaccharideResidue.from_monosaccharide(base)
            self.assertAlmostEqual(base.mass() - water_mass, residue.mass(), 3)

    def test_clone_coherence(self):
        for base in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[base]
            residue = glycan_composition.FrozenMonosaccharideResidue.from_monosaccharide(base)
            dup = residue.clone()
            self.assertAlmostEqual(dup.mass(), residue.mass(), 3)

    def test_derivatize(self):
        reference = {"GlcNAc": 245.1263, "Fuc": 174.0892, "NeuAc": 361.1737}
        with self.assertRaises(glycan_composition.FrozenError):
            for name in ["GlcNAc", "NeuAc", "Fuc"]:
                base = monosaccharides[name]
                residue = glycan_composition.FrozenMonosaccharideResidue.from_monosaccharide(base)
                composition_transform.derivatize(base, "methyl")
                composition_transform.derivatize(residue, "methyl")
                self.assertAlmostEqual(
                    base.mass() - (water_mass + 2 * (Substituent("methyl").mass() - glypy.Composition("H").mass * 2)),
                    residue.mass(), 3)
                self.assertAlmostEqual(reference[name], residue.mass(), 3)

                base = monosaccharides[name]
                composition_transform.derivatize(base, "methyl")
                residue = glycan_composition.FrozenMonosaccharideResidue.from_monosaccharide(base)

                self.assertAlmostEqual(
                    base.mass() - (water_mass + 2 * (Substituent("methyl").mass() - glypy.Composition("H").mass * 2)),
                    residue.mass(), 3)

                self.assertAlmostEqual(reference[name], residue.mass(), 3)

    def test_acidic_special_case(self):
        acidic_hexose1 = glycan_composition.FrozenMonosaccharideResidue.from_iupac_lite("HexA")
        acidic_hexose2 = glycan_composition.FrozenMonosaccharideResidue.from_iupac_lite("aHex")
        hexose = glycan_composition.FrozenMonosaccharideResidue.from_iupac_lite("Hex")
        self.assertAlmostEqual(acidic_hexose1.mass() - hexose.mass(), 13.9792, 3)
        self.assertAlmostEqual(acidic_hexose2.mass() - hexose.mass(), 13.9792, 3)


class GlycanCompositionTests(unittest.TestCase):
    GlycanCompositionType = GlycanComposition

    def test_from_glycan(self):
        glyc = glycans["N-Linked Core"]
        comp = self.GlycanCompositionType.from_glycan(glyc)
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)

    def test_derivatize(self):
        glyc = glycans["N-Linked Core"]
        composition_transform.derivatize(glyc, "methyl")
        comp = self.GlycanCompositionType.from_glycan(glyc)
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)

        comp = self.GlycanCompositionType.from_glycan(glycans["N-Linked Core"])
        composition_transform.derivatize(comp, "methyl")
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)
        self.assertAlmostEqual(self.GlycanCompositionType.parse(comp).mass(), comp.mass(), 3)

    def test_contains(self):
        glyc = glycans["N-Linked Core"]
        comp = self.GlycanCompositionType.from_glycan(glyc)
        self.assertTrue("Man" in comp)

    def test_serialize(self):
        ref = '{Man:3; Glc2NAc:2}'
        glyc = glycans["N-Linked Core"]
        comp = self.GlycanCompositionType.from_glycan(glyc)
        self.assertEqual(comp.serialize(), ref)

        ref2 = '{Hex:3; Hex2NAc:2}'
        comp.drop_stems()
        self.assertEqual(comp.serialize(), ref2)

    def test_parse(self):
        ref = '{Man:3; Glc2NAc:2}'
        glyc = glycans["N-Linked Core"]
        comp = self.GlycanCompositionType.from_glycan(glyc)

        comp2 = self.GlycanCompositionType.parse(ref)
        self.assertEqual(comp, comp2)

    def test_update(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = self.GlycanCompositionType.parse(ref)
        self.assertEqual(comp["Man"], 3)
        comp.update({"Man": 1})
        self.assertEqual(comp["Man"], 1)
        comp.update(Man=2)
        self.assertEqual(comp["Man"], 2)

    def test_arithmetic(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = self.GlycanCompositionType.parse(ref)

        comp["Man"] += 5
        self.assertEqual(comp["Man"], 8)
        comp2 = self.GlycanCompositionType(Neu5NAc=2, Fuc=1, Man=1)
        self.assertEqual(comp2 + comp, self.GlycanCompositionType(
            Neu5NAc=2, Fuc=1, Man=9, Glc2NAc=2))
        comp -= {'Man': 5}
        self.assertEqual(comp, self.GlycanCompositionType(
            Man=4, Glc2NAc=2) - self.GlycanCompositionType(Man=1, Glc2NAc=0))
        comp *= 2
        self.assertEqual(comp, (self.GlycanCompositionType(
            Man=4, Glc2NAc=2) - self.GlycanCompositionType(Man=1, Glc2NAc=0)) * 2)

        comp3 = comp2.clone()
        comp3 += comp2
        self.assertEqual(comp3, comp2 * 2)

    def test_total_composition(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = self.GlycanCompositionType.parse(ref)
        glyc = glycans["N-Linked Core"]
        self.assertEqual(comp.total_composition(), glyc.total_composition())

    def test_drops(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = self.GlycanCompositionType.parse(ref)
        self.assertEqual(comp["Man"], 3)
        self.assertEqual(comp["Glc2NAc"], 2)
        comp.drop_positions()
        self.assertEqual(comp["GlcNAc"], 2)
        self.assertEqual(comp["Glc2NAc"], 0)

        comp.drop_stems()
        self.assertEqual(comp["HexNAc"], 2)
        self.assertEqual(comp["GlcNAc"], 0)
        self.assertEqual(comp["Glc2NAc"], 0)
        self.assertEqual(comp["Hex"], 3)
        self.assertEqual(comp["Man"], 0)

    def test_parse_derivatized_reduced(self):
        x = self.GlycanCompositionType.parse('{Fuc^Me:1; Hex^Me:5; HexNAc^Me:4; Neu5NAc^Me:1}$C1H4')
        self.assertAlmostEqual(x.mass(), 2598.3402, 4)


class FrozenGlycanCompositionTests(GlycanCompositionTests):
    GlycanCompositionType = glycan_composition.FrozenGlycanComposition

    def test_derivatize(self):
        glyc = glycans["N-Linked Core"]
        with self.assertRaises(glycan_composition.FrozenError):
            composition_transform.derivatize(glyc, "methyl")
            comp = self.GlycanCompositionType.from_glycan(glyc)
            self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)

            comp = self.GlycanCompositionType.from_glycan(glycans["N-Linked Core"])
            composition_transform.derivatize(comp, "methyl")
            self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)
            self.assertAlmostEqual(self.GlycanCompositionType.parse(comp).mass(), comp.mass(), 3)

    def test_parse_derivatized(self):
        glyc = glycans["N-Linked Core"]
        composition_transform.derivatize(glyc, "methyl")
        comp = GlycanComposition.from_glycan(glycans["N-Linked Core"])
        composition_transform.derivatize(comp, "methyl")
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)
        frozen_comp = self.GlycanCompositionType.parse(comp)
        self.assertAlmostEqual(glyc.mass(), frozen_comp.mass(), 3)

    def test_serialize(self):
        ref = '{Man:3; Glc2NAc:2}'
        glyc = glycans["N-Linked Core"]
        comp = self.GlycanCompositionType.from_glycan(glyc)
        self.assertEqual(comp.serialize(), ref)

    def test_drops(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = self.GlycanCompositionType.parse(ref)
        self.assertEqual(comp["Man"], 3)
        self.assertEqual(comp["Glc2NAc"], 2)
        with self.assertRaises(glycan_composition.FrozenError):
            comp.drop_positions()


class SubstituentResidueTests(unittest.TestCase):
    def test_parse(self):
        residue = glycan_composition.from_iupac_lite("@n_acetyl")
        hexose = MonosaccharideResidue.from_iupac_lite("Hex")
        hexnac = MonosaccharideResidue.from_iupac_lite("HexNAc")
        self.assertAlmostEqual(residue.mass(), hexnac.mass() - hexose.mass(), 3)

    def test_to_iupac_lite(self):
        reference = "@n_acetyl"
        residue = glycan_composition.from_iupac_lite(reference)
        self.assertEqual(reference, str(residue))

    def test_hashable(self):
        c = dict()
        residue = glycan_composition.from_iupac_lite("@n_acetyl")
        c[residue] = 5
        self.assertEqual(c[glycan_composition.from_iupac_lite("@n_acetyl")], 5)

    def test_equality(self):
        residue = glycan_composition.from_iupac_lite("@n_acetyl")
        self.assertEqual(residue, glycan_composition.from_iupac_lite("@n_acetyl"))
        self.assertNotEqual(residue, "n_acetyl")


if __name__ == '__main__':
    unittest.main()
