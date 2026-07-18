import json,math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import ValidationError
from mesim.thermo.ideal import load_correlations
from mesim.thermo.steam_tables import load_steam_tables_data,steam_saturation_pressure_pa,steam_tables_fugacity_coefficients,steam_tables_tp_flash
from mesim.thermo.systems import STEAM_TABLES_WATER,SteamTablesSystem,create_thermo_system
class SteamTablesTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.path=ROOT/"data/correlations/steam-tables-v1.json";cls.data=load_steam_tables_data(cls.path);cls.water=next(x for x in load_correlations(ROOT/"data/correlations/ideal-v1.json") if x.compound_id=="Water")
 def test_region4_source_is_frozen(self):
  self.assertEqual(len(self.data.region4_coefficients),10);self.assertEqual(self.data.iapws_source_sha256,"b5ed91b55175783df46df0670b4cb607df5cb80170ad74ef07d744d8efb9837f")
  self.assertTrue(math.isclose(steam_saturation_pressure_pa(self.data,450),932041.0791359359,rel_tol=1e-15))
 def test_repeatable_direct_class_golden_parity(self):
  g=json.loads((ROOT/"tests/golden/steam-tables-water-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/steam-tables-water-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r)
  self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.SteamTablesPropertyPackage");self.assertEqual(g["source"]["property_package_construction"],"direct-class-over-case-compound-domain")
  i=g["inputs"];o=g["outputs"];self.assertEqual(steam_tables_fugacity_coefficients(self.data,self.water,(1,),i["temperature_k"],i["pressure_pa"],"liquid"),tuple(o["liquid_fugacity_coefficients"]));self.assertEqual(steam_tables_fugacity_coefficients(self.data,self.water,(1,),i["temperature_k"],i["pressure_pa"],"vapor"),tuple(o["vapor_fugacity_coefficients"]))
  result=steam_tables_tp_flash(self.data,self.water,(1,),i["temperature_k"],i["pressure_pa"]);self.assertEqual((result.liquid_fraction,result.vapor_fraction,result.liquid_composition,result.vapor_composition,result.iterations),(o["liquid_fraction"],o["vapor_fraction"],tuple(o["liquid_composition"]),tuple(o["vapor_composition"]),o["iterations"]))
 def test_system_and_rejections(self):
  system=create_thermo_system(STEAM_TABLES_WATER,data=self.data,water=self.water);self.assertIsInstance(system,SteamTablesSystem);self.assertEqual(system.tp_flash((1,),450,1e6),steam_tables_tp_flash(self.data,self.water,(1,),450,1e6))
  with self.assertRaises(ValidationError):steam_saturation_pressure_pa(self.data,700)
  with self.assertRaises(ValidationError):steam_tables_fugacity_coefficients(self.data,self.water,(.5,.5),450,1e6,"liquid")
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["region4_coefficients"].pop()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_steam_tables_data(p)
if __name__=="__main__":unittest.main()
