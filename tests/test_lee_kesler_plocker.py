import json,math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import MissingCompoundData,ValidationError
from mesim.thermo.lee_kesler_plocker import load_lkp_data,lkp_fugacity_coefficients,lkp_tp_flash
from mesim.thermo.systems import LEE_KESLER_PLOCKER_METHANE_N_PENTANE,LeeKeslerPlockerSystem,create_thermo_system
class LeeKeslerPlockerTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.path=ROOT/"data/interactions/lkp-v1.json";cls.data=load_lkp_data(cls.path);cls.ids=("Methane","N-pentane")
 def test_complete_interaction_domain_and_scoped_basis_are_frozen(self):
  self.assertEqual(len(self.data.interaction_pairs),140);self.assertEqual(self.data.interaction_resource_sha256,"20e48a501f4d568a1a41b6d265c3d5d7bdaddc729e9aeedce1ce4e3961a379ed")
  self.assertEqual(self.data.interaction("Methane","N-pentane"),1.2401);self.assertEqual(self.data.interaction("N-pentane","Methane"),1.2401);self.assertEqual(self.data.interaction("Methane","Methane"),1.0)
  self.assertEqual(tuple(x.critical_volume for x in self.data.compounds),(0.0986,0.311))
 def test_repeatable_phase_fugacity_parity(self):
  g=json.loads((ROOT/"tests/golden/lkp-methane-n-pentane-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/lkp-methane-n-pentane-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r);self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.LKPPropertyPackage")
  i=g["inputs"];o=g["outputs"]
  for phase,key in (("liquid","liquid_fugacity_coefficients"),("vapor","vapor_fugacity_coefficients")):
   actual=lkp_fugacity_coefficients(self.data,self.ids,tuple(i["composition"]),i["temperature_k"],i["pressure_pa"],phase)
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,abs_tol=1e-9))
 def test_tighter_flash_matches_saved_reference(self):
  g=json.loads((ROOT/"tests/golden/lkp-methane-n-pentane-state.json").read_text(encoding="utf-8-sig"));i=g["inputs"];o=g["outputs"];result=lkp_tp_flash(self.data,self.ids,tuple(i["composition"]),i["temperature_k"],i["pressure_pa"])
  self.assertTrue(math.isclose(result.liquid_fraction,o["liquid_fraction"],abs_tol=1e-4));self.assertTrue(math.isclose(result.vapor_fraction,o["vapor_fraction"],abs_tol=1e-4))
  for key,actual in (("liquid_composition",result.liquid_composition),("vapor_composition",result.vapor_composition),("equilibrium_ratios",result.equilibrium_ratios)):
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=1e-4,abs_tol=1e-4))
 def test_system_and_invalid_data(self):
  system=create_thermo_system(LEE_KESLER_PLOCKER_METHANE_N_PENTANE,data=self.data,compound_ids=self.ids);self.assertIsInstance(system,LeeKeslerPlockerSystem);self.assertEqual(system.tp_flash((.5,.5),350,1e6),lkp_tp_flash(self.data,self.ids,(.5,.5),350,1e6))
  with self.assertRaises(MissingCompoundData):self.data.compound("methane")
  with self.assertRaises(ValidationError):lkp_fugacity_coefficients(self.data,self.ids,(.5,.6),350,1e6,"vapor")
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["interaction_pairs"].pop()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_lkp_data(p)
if __name__=="__main__":unittest.main()
