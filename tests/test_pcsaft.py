import json,math,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import ValidationError
from mesim.thermo.pcsaft import load_pcsaft_data,pcsaft_fugacity_coefficients,pcsaft_tp_flash
from mesim.thermo.systems import PCSAFT_METHANE_ETHANE,PCSAFTSystem,create_thermo_system
class PCSAFTTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.path=ROOT/"data/interactions/pcsaft-v1.json";cls.data=load_pcsaft_data(cls.path)
 def test_complete_source_tables_and_scoped_parameters(self):
  self.assertEqual(len(self.data.pure_parameters),94);self.assertEqual(len(self.data.interaction_parameters),33);self.assertEqual(len({x.cas_number for x in self.data.pure_parameters}),93);self.assertEqual(self.data.eos_source_sha256,"c28c9aff3310c331927d50ff6e02825ed0c26efe7c6fd64f0fc5e7c82598508e")
  methane=self.data.parameter("74-82-8");ethane=self.data.parameter("74-84-0");self.assertEqual((methane.segment_count,methane.segment_diameter_angstrom,methane.dispersion_energy_k),(1,3.7039,150.03));self.assertEqual((ethane.segment_count,ethane.segment_diameter_angstrom,ethane.dispersion_energy_k),(1.6069,3.5206,191.42));self.assertEqual(self.data.kij("74-82-8","74-84-0"),0)
 def test_repeatable_golden_equation_and_flash_parity(self):
  g=json.loads((ROOT/"tests/golden/pc-saft-methane-ethane-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/pc-saft-methane-ethane-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r);self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.AdvancedEOS.PCSAFT2PropertyPackage")
  i=g["inputs"];o=g["outputs"]
  for phase,key in (("liquid","liquid_fugacity_coefficients"),("vapor","vapor_fugacity_coefficients")):
   actual=pcsaft_fugacity_coefficients(self.data,i["compounds"],i["composition"],i["temperature_k"],i["pressure_pa"],phase)
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=2e-10))
  x=tuple(o["liquid_composition"]);y=tuple(o["vapor_composition"]);phil=pcsaft_fugacity_coefficients(self.data,i["compounds"],x,i["temperature_k"],i["pressure_pa"],"liquid");phiv=pcsaft_fugacity_coefficients(self.data,i["compounds"],y,i["temperature_k"],i["pressure_pa"],"vapor");reference_residual=max(abs(math.log(x[j]*phil[j]/(y[j]*phiv[j]))) for j in range(2));self.assertLess(reference_residual,7e-5)
  result=pcsaft_tp_flash(self.data,i["compounds"],i["composition"],i["temperature_k"],i["pressure_pa"]);self.assertLess(result.residual_norm,1e-7)
  for value,expected in ((result.liquid_fraction,o["liquid_fraction"]),(result.vapor_fraction,o["vapor_fraction"])):self.assertTrue(math.isclose(value,expected,rel_tol=0,abs_tol=7e-5))
  for actual,key in ((result.liquid_composition,"liquid_composition"),(result.vapor_composition,"vapor_composition"),(result.equilibrium_ratios,"equilibrium_ratios")):
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=0,abs_tol=7e-5))
 def test_system_and_rejections(self):
  system=create_thermo_system(PCSAFT_METHANE_ETHANE,data=self.data);self.assertIsInstance(system,PCSAFTSystem);self.assertEqual(system.tp_flash((.7,.3),220,3e6),pcsaft_tp_flash(self.data,self.data.compound_ids,(.7,.3),220,3e6))
  with self.assertRaises(ValidationError):pcsaft_fugacity_coefficients(self.data,tuple(reversed(self.data.compound_ids)),(.7,.3),220,3e6,"vapor")
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["pure_parameters"].pop()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_pcsaft_data(p)
if __name__=="__main__":unittest.main()
