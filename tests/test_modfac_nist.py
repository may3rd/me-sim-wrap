import json, math, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import MissingCompoundData, ValidationError
from mesim.thermo.modfac import load_modfac_data, modfac_activity_coefficients
from mesim.thermo.systems import MODFAC_NIST_1_PROPANOL_WATER, ModfacNistSystem, create_thermo_system

class ModfacNistTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.path=ROOT/"data/interactions/modfac-nist-v1.json";cls.data=load_modfac_data(cls.path)
 def test_full_nist_domains_are_frozen(self):
  self.assertEqual(self.data.model,"Modified UNIFAC (NIST)");self.assertEqual(len(self.data.groups),201);self.assertEqual(len(self.data.interaction_pairs),1969)
  self.assertEqual(self.data.groups_sha256,"46d652018b8205cab3274d237cc0a515127773a8a684082fa46f825ab60457ff")
  self.assertEqual(self.data.interactions_sha256,"6074a15684c188e48a8783f44194462dd91200fd65c6a7d6e0c978bbea8d4f7d")
  p=self.data.compound("1-propanol");self.assertEqual((p.q,p.r),(3.3697,3.1277));self.assertEqual(self.data.coefficients(1,7),(1391.3,-3.6156,0.0011439999999999998));self.assertEqual(self.data.coefficients(7,1),(-17.25,0.8389,0.0009021))
 def test_repeatable_dwsim_activity_parity(self):
  g=json.loads((ROOT/"tests/golden/modfac-nist-1-propanol-water-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/modfac-nist-1-propanol-water-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r)
  self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.NISTMFACPropertyPackage")
  i=g["inputs"];a=modfac_activity_coefficients(self.data,tuple(i["compounds"]),tuple(i["composition"]),i["temperature_k"])
  self.assertEqual(a,tuple(g["outputs"]["activity_coefficients"]))
 def test_registered_system(self):
  s=create_thermo_system(MODFAC_NIST_1_PROPANOL_WATER,data=self.data,compound_ids=("1-propanol","Water"));self.assertIsInstance(s,ModfacNistSystem);self.assertEqual(s.activity_coefficients((.5,.5),350),modfac_activity_coefficients(self.data,s.compound_ids,(.5,.5),350))
  with self.assertRaises(ValidationError):create_thermo_system(MODFAC_NIST_1_PROPANOL_WATER,data=load_modfac_data(ROOT/"data/interactions/modfac-dortmund-v1.json"),compound_ids=("1-propanol","Water"))
 def test_invalid_inputs_rejected(self):
  with self.assertRaises(MissingCompoundData):self.data.compound("1-Propanol")
  with self.assertRaises(ValidationError):modfac_activity_coefficients(self.data,("1-propanol","Water"),(.5,.6),350)
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["interaction_pairs"].pop()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_modfac_data(p)
if __name__=="__main__":unittest.main()
