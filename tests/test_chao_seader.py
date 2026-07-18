import json, math, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from mesim.errors import MissingCompoundData, ValidationError
from mesim.thermo.chao_seader import load_chao_seader_data, chao_seader_liquid_fugacity_coefficients, chao_seader_vapor_fugacity_coefficients, chao_seader_tp_flash
from mesim.thermo.systems import CHAO_SEADER_METHANE_N_PENTANE, ChaoSeaderSystem, create_thermo_system

class ChaoSeaderTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.path=ROOT/"data/interactions/chao-seader-v1.json";cls.data=load_chao_seader_data(cls.path);cls.ids=("Methane","N-pentane")
 def test_scoped_source_data_are_frozen(self):
  self.assertEqual(self.data.model,"Chao-Seader")
  self.assertEqual(tuple(x.compound_id for x in self.data.compounds),self.ids)
  self.assertEqual(self.data.runtime_assembly_sha256,"c5a038c86cdfd3304be5e283527ab377df7579ca306eac69950ab4b8c495c544")
  self.assertEqual(self.data.property_package_source_sha256,"284ae761082d3fab4104d0c636a1023349cee15719757c3d932a7118478f023f")
  self.assertEqual(self.data.model_source_sha256,"5fe15931118e56a149d7b321fa6120e94b37aea1c763f541b54f760c2640d405")
  self.assertEqual((self.data.compound("Methane").liquid_molar_volume,self.data.compound("N-pentane").solubility_parameter),(37.8392,0.00342968764086))
 def test_repeatable_phase_fugacity_parity(self):
  g=json.loads((ROOT/"tests/golden/chao-seader-methane-n-pentane-state.json").read_text(encoding="utf-8-sig"));r=json.loads((ROOT/"tests/golden/chao-seader-methane-n-pentane-state-repeat.json").read_text(encoding="utf-8-sig"));self.assertEqual(g,r)
  self.assertEqual(g["source"]["property_package_class"],"DWSIM.Thermodynamics.PropertyPackages.ChaoSeaderPropertyPackage")
  i=g["inputs"];o=g["outputs"]
  liquid=chao_seader_liquid_fugacity_coefficients(self.data,self.ids,tuple(i["composition"]),i["temperature_k"],i["pressure_pa"]);vapor=chao_seader_vapor_fugacity_coefficients(self.data,self.ids,tuple(i["composition"]),i["temperature_k"],i["pressure_pa"])
  for actual,expected in zip(liquid,o["liquid_fugacity_coefficients"]):self.assertTrue(math.isclose(actual,expected,rel_tol=3e-15))
  for actual,expected in zip(vapor,o["vapor_fugacity_coefficients"]):self.assertTrue(math.isclose(actual,expected,rel_tol=3e-15))
 def test_tighter_flash_matches_reference_with_documented_tolerance(self):
  g=json.loads((ROOT/"tests/golden/chao-seader-methane-n-pentane-state.json").read_text(encoding="utf-8-sig"));i=g["inputs"];o=g["outputs"];result=chao_seader_tp_flash(self.data,self.ids,tuple(i["composition"]),i["temperature_k"],i["pressure_pa"])
  self.assertTrue(math.isclose(result.liquid_fraction,o["liquid_fraction"],abs_tol=1e-4));self.assertTrue(math.isclose(result.vapor_fraction,o["vapor_fraction"],abs_tol=1e-4))
  for key,actual in (("liquid_composition",result.liquid_composition),("vapor_composition",result.vapor_composition),("equilibrium_ratios",result.equilibrium_ratios)):
   for value,expected in zip(actual,o[key]):self.assertTrue(math.isclose(value,expected,rel_tol=1e-4,abs_tol=1e-4))
  self.assertLess(result.iterations,100)
 def test_registered_system_and_rejections(self):
  system=create_thermo_system(CHAO_SEADER_METHANE_N_PENTANE,data=self.data,compound_ids=self.ids);self.assertIsInstance(system,ChaoSeaderSystem);self.assertEqual(system.tp_flash((.5,.5),350,1e6),chao_seader_tp_flash(self.data,self.ids,(.5,.5),350,1e6))
  with self.assertRaises(MissingCompoundData):self.data.compound("methane")
  with self.assertRaises(ValidationError):chao_seader_liquid_fugacity_coefficients(self.data,self.ids,(.5,.6),350,1e6)
  d=json.loads(self.path.read_text(encoding="utf-8-sig"));d["compounds"].pop()
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/"bad.json";p.write_text(json.dumps(d),encoding="utf-8")
   with self.assertRaises(ValidationError):load_chao_seader_data(p)
if __name__=="__main__":unittest.main()
