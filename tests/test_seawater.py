import json,math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import ValidationError
from mesim.thermo.seawater import load_seawater_data,seawater_fugacity_coefficients,seawater_salinity,seawater_tp_flash,seawater_vapor_pressure_pa
from mesim.thermo.systems import SEAWATER_WATER_SALT,SeawaterSystem,create_thermo_system

class SeawaterTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.path=ROOT/"data/correlations/seawater-v1.json";cls.data=load_seawater_data(cls.path)
 def test_source_and_two_salinity_paths_are_frozen(self):
  self.assertEqual(self.data.compound_ids,("Water","Salt"));self.assertEqual(self.data.molecular_weights,(18.01528,58.44));self.assertEqual(len(self.data.vapor_pressure_coefficients),6)
  self.assertEqual(self.data.property_package_source_sha256,"435448de7115d5c00974f80e284ed7cef5976491a0f4974c7502a7c548739f79")
  self.assertTrue(math.isclose(seawater_salinity(self.data,(.98,.02)),self.data.probe_salinity,rel_tol=0,abs_tol=1e-16));self.assertEqual(self.data.current_stream_salinity,.12)
  self.assertEqual(seawater_vapor_pressure_pa(self.data,350,composition=(.98,.02)),self.data.probe_vapor_pressure_pa)
 def test_repeatable_direct_class_golden_parity(self):
  g=json.loads((ROOT/"tests/golden/seawater-water-salt-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/seawater-water-salt-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r)
  self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.SeawaterPropertyPackage");self.assertEqual(g["source"]["property_package_construction"],"direct-class-over-case-compound-domain")
  i=g["inputs"];o=g["outputs"]
  for phase,key in (("liquid","liquid_fugacity_coefficients"),("vapor","vapor_fugacity_coefficients")):
   self.assertEqual(seawater_fugacity_coefficients(self.data,i["composition"],i["temperature_k"],i["pressure_pa"],phase),tuple(o[key]))
  result=seawater_tp_flash(self.data,i["composition"],i["temperature_k"],i["pressure_pa"])
  self.assertEqual((result.liquid_fraction,result.vapor_fraction,result.liquid_composition,result.vapor_composition,result.iterations,result.equilibrium_ratios),(o["liquid_fraction"],o["vapor_fraction"],tuple(o["liquid_composition"]),tuple(o["vapor_composition"]),o["iterations"],tuple(o["equilibrium_ratios"])))
 def test_system_and_rejections(self):
  system=create_thermo_system(SEAWATER_WATER_SALT,data=self.data);self.assertIsInstance(system,SeawaterSystem);self.assertEqual(system.tp_flash((.98,.02),350,1e5),seawater_tp_flash(self.data,(.98,.02),350,1e5))
  with self.assertRaises(ValidationError):seawater_salinity(self.data,(1,))
  with self.assertRaises(ValidationError):seawater_vapor_pressure_pa(self.data,350)
  with self.assertRaises(ValidationError):seawater_tp_flash(self.data,(.98,.02),350,1000)
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["compound_ids"].reverse()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_seawater_data(p)
if __name__=="__main__":unittest.main()
