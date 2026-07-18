import json,math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import ValidationError
from mesim.thermo.ideal import load_correlations
from mesim.thermo.ideal_electrolyte import load_ideal_electrolyte_data,ideal_electrolyte_fugacity_coefficients,ideal_electrolyte_molalities,ideal_electrolyte_tp_flash
from mesim.thermo.systems import IDEAL_AQUEOUS_ELECTROLYTE_WATER_SODIUM_CHLORIDE,IdealAqueousElectrolyteSystem,create_thermo_system
class IdealElectrolyteTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.path=ROOT/"data/correlations/ideal-aqueous-electrolyte-v1.json";cls.data=load_ideal_electrolyte_data(cls.path);cls.water=next(x for x in load_correlations(ROOT/"data/correlations/ideal-v1.json") if x.compound_id=="Water")
 def test_source_domain_and_molality_are_frozen(self):
  self.assertEqual(self.data.compound_ids,("Water","Sodium (ion)","Chloride (ion)"));self.assertEqual(tuple(x.formula for x in self.data.compounds),("HOH","Na+","Cl-"));self.assertEqual(self.data.database_source_sha256,"bc72f2a5e4ac273a2ea4d881b4fad2bc5d0e7e9f470f019751205a00004c7fb1")
  actual=ideal_electrolyte_molalities(self.data,self.data.compound_ids,self.data.probe_composition)
  for value,expected in zip(actual,self.data.probe_molalities):self.assertTrue(math.isclose(value,expected,rel_tol=2e-15))
 def test_repeatable_saved_case_golden_parity(self):
  g=json.loads((ROOT/"tests/golden/ideal-aqueous-electrolyte-water-sodium-chloride-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/ideal-aqueous-electrolyte-water-sodium-chloride-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r);self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.IdealElectrolytePropertyPackage")
  i=g["inputs"];o=g["outputs"]
  for phase,key in (("liquid","liquid_fugacity_coefficients"),("vapor","vapor_fugacity_coefficients")):
   actual=ideal_electrolyte_fugacity_coefficients(self.data,self.water,i["compounds"],i["composition"],i["temperature_k"],i["pressure_pa"],phase)
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=2e-15))
  result=ideal_electrolyte_tp_flash(self.data,self.water,i["compounds"],i["composition"],i["temperature_k"],i["pressure_pa"]);self.assertEqual((result.liquid_fraction,result.vapor_fraction,result.liquid_composition,result.vapor_composition,result.equilibrium_ratios,result.iterations),(o["liquid_fraction"],o["vapor_fraction"],tuple(o["liquid_composition"]),tuple(o["vapor_composition"]),tuple(o["equilibrium_ratios"]),o["iterations"]))
 def test_system_and_rejections(self):
  system=create_thermo_system(IDEAL_AQUEOUS_ELECTROLYTE_WATER_SODIUM_CHLORIDE,data=self.data,water=self.water);self.assertIsInstance(system,IdealAqueousElectrolyteSystem);self.assertEqual(system.tp_flash((.98,.01,.01),298.15,1e5),ideal_electrolyte_tp_flash(self.data,self.water,self.data.compound_ids,(.98,.01,.01),298.15,1e5))
  with self.assertRaises(ValidationError):ideal_electrolyte_molalities(self.data,self.data.compound_ids,(.98,.02,0))
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["compounds"][1]["is_ion"]=False
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_ideal_electrolyte_data(p)
if __name__=="__main__":unittest.main()
