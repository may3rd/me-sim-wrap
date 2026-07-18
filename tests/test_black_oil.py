import json,math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import ValidationError
from mesim.thermo.black_oil import load_black_oil_data,black_oil_component_vaporized_fraction,black_oil_fugacity_coefficients,black_oil_tp_flash,black_oil_vapor_pressure_pa
from mesim.thermo.steam_tables import load_steam_tables_data
from mesim.thermo.systems import BLACK_OIL_N_PENTANE_N_HEXANE,BlackOilSystem,create_thermo_system

class BlackOilTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.path=ROOT/"data/correlations/black-oil-v1.json";cls.data=load_black_oil_data(cls.path);cls.steam=load_steam_tables_data(ROOT/"data/correlations/steam-tables-v1.json")
 def test_source_and_runtime_probe_are_frozen(self):
  self.assertEqual(self.data.compound_ids,("n-Pentane","n-Hexane"));self.assertEqual(self.data.case_sha256,"e734c1178212f7c7eac1404c61cfabc7a11a226a11976548ba03d9b30b0b2cec");self.assertEqual(self.data.flash_source_sha256,"ab6136e745474fc0ee702f8f7848fcc6784430c96d633752c0015a50e61c06ee")
  pressures=tuple(black_oil_vapor_pressure_pa(self.data,self.steam,name,self.data.probe_temperature_k) for name in self.data.compound_ids);fractions=tuple(black_oil_component_vaporized_fraction(record,self.data.probe_temperature_k,self.data.probe_pressure_pa) for record in self.data.compounds)
  for actual,expected in zip(pressures,self.data.probe_vapor_pressures_pa):self.assertTrue(math.isclose(actual,expected,rel_tol=2e-15))
  for actual,expected in zip(fractions,self.data.probe_vaporized_fractions):self.assertTrue(math.isclose(actual,expected,rel_tol=2e-15))
 def test_repeatable_saved_case_golden_parity(self):
  g=json.loads((ROOT/"tests/golden/black-oil-n-pentane-n-hexane-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/black-oil-n-pentane-n-hexane-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r);self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.BlackOilPropertyPackage");self.assertEqual(g["source"]["property_package_construction"],"deserialized-from-case")
  i=g["inputs"];o=g["outputs"]
  for phase,key in (("liquid","liquid_fugacity_coefficients"),("vapor","vapor_fugacity_coefficients")):
   actual=black_oil_fugacity_coefficients(self.data,self.steam,i["compounds"],i["composition"],i["temperature_k"],i["pressure_pa"],phase)
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=2e-15,abs_tol=1e-20))
  result=black_oil_tp_flash(self.data,i["compounds"],i["composition"],i["temperature_k"],i["pressure_pa"])
  for value,expected in ((result.liquid_fraction,o["liquid_fraction"]),(result.vapor_fraction,o["vapor_fraction"])):self.assertTrue(math.isclose(value,expected,rel_tol=2e-15))
  for actual,key in ((result.liquid_composition,"liquid_composition"),(result.vapor_composition,"vapor_composition"),(result.equilibrium_ratios,"equilibrium_ratios")):
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=2e-15,abs_tol=1e-15))
  self.assertEqual(result.iterations,o["iterations"])
 def test_system_and_rejections(self):
  system=create_thermo_system(BLACK_OIL_N_PENTANE_N_HEXANE,data=self.data,steam_data=self.steam);self.assertIsInstance(system,BlackOilSystem);self.assertEqual(system.tp_flash((.6,.4),350,1e6),black_oil_tp_flash(self.data,self.data.compound_ids,(.6,.4),350,1e6))
  with self.assertRaises(ValidationError):black_oil_tp_flash(self.data,tuple(reversed(self.data.compound_ids)),(.6,.4),350,1e6)
  with self.assertRaises(ValidationError):black_oil_vapor_pressure_pa(self.data,self.steam,"Methane",350)
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["compounds"][0]["specific_gravity_oil"]=1.2
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_black_oil_data(p)
if __name__=="__main__":unittest.main()
