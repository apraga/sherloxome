{ pkgs, lib, python3, fetchPypi }:

let
  my-python-packages = p: with p; [
    biopython
    matplotlib
    pkginfo
    numpy
    pyyaml
    scipy
    pysam
    frozendict
    pip # declared as a runtime dependency in NEAT's pyproject.toml
  ];
in
python3.pkgs.buildPythonApplication rec {
  pname = "neat-genreads";
  version = "4.7.0";
  pyproject = true;

  src = fetchPypi {
    inherit version;
    pname = "neat_genreads";
    hash = "sha256-WqI+2CxcYkEPfqsdhSZMN2mwnd4RikBfLes/Lnc/9wU=";
  };

  build-system = [ python3.pkgs.poetry-core ];

  dependencies = my-python-packages python3.pkgs;

  # NEAT's pyproject.toml pins pip >=25.2, newer than the nixpkgs-provided
  # pip; pip isn't actually used at runtime by NEAT's code, so relax it.
  pythonRelaxDeps = [ "pip" ];

  # No test suite is shipped in the PyPI sdist.
  doCheck = false;

  pythonImportsCheck = [ "neat" ];

  meta = with lib; {
    description = "NGS read/variant simulation toolkit (NEAT)";
    homepage = "https://github.com/ncsa/NEAT";
    license = licenses.bsd3;
    mainProgram = "neat";
  };
}
