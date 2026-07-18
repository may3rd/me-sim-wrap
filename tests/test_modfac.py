import json, math, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import MissingCompoundData, ValidationError
from mesim.thermo.modfac import load_modfac_data, modfac_activity_coefficients
from mesim.thermo.systems import MODFAC_DORTMUND_1_PROPANOL_WATER, ModfacDortmundSystem, create_thermo_system

class ModfacTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.path=ROOT/"data/interactions/modfac-dortmund-v1.json";cls.data=load_modfac_data(cls.path)
 def test_full_dortmund_domains_are_frozen(self):
  self.assertEqual(len(self.data.groups),108);self.assertEqual(len(self.data.interaction_pairs),1167)
  self.assertEqual(self.data.groups_sha256,"275b73f9f6524156b1a733d1d921e9af065ada901c4dc387d90387b95f7fac95")
  self.assertEqual(self.data.interactions_sha256,"bc50dd0892752ad5214534b627b29c1fbb8fd875d8634002b6d02d7cce744439")
  p=self.data.compound("1-propanol");self.assertEqual((p.q,p.r),(3.3697,3.1277));self.assertEqual(self.data.coefficients(1,7),(1391.3,-3.6156,0.001144));self.assertEqual(self.data.coefficients(7,1),(-17.253,0.8389,0.0009021))
 def test_repeatable_dwsim_activity_parity(self):
  g=json.loads((ROOT/"tests/golden/modfac-dortmund-1-propanol-water-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/modfac-dortmund-1-propanol-water-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r)
  self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.MODFACPropertyPackage")
  i=g["inputs"];a=modfac_activity_coefficients(self.data,tuple(i["compounds"]),tuple(i["composition"]),i["temperature_k"])
  for x,y in zip(a,g["outputs"]["activity_coefficients"]):self.assertTrue(math.isclose(x,y,rel_tol=3e-15))
 def test_registered_system(self):
  s=create_thermo_system(MODFAC_DORTMUND_1_PROPANOL_WATER,data=self.data,compound_ids=("1-propanol","Water"));self.assertIsInstance(s,ModfacDortmundSystem);self.assertEqual(s.activity_coefficients((.5,.5),350),modfac_activity_coefficients(self.data,s.compound_ids,(.5,.5),350))
 def test_invalid_inputs_rejected(self):
  with self.assertRaises(MissingCompoundData):self.data.compound("1-Propanol")
  with self.assertRaises(ValidationError):modfac_activity_coefficients(self.data,("1-propanol","Water"),(.5,.6),350)
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["groups"].pop()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_modfac_data(p)
if __name__=="__main__":unittest.main()
