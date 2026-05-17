{ config, lib, pkgs, fetchFromGitHub, stdenv  }:

stdenv.mkDerivation {
  pname = "simuscop";
  version = "0.1.0";

  src = fetchFromGitHub {
    owner = "qasimyu";
    repo = "simuscop";
    rev = "590bff10e640bca3b520e4737aad9bc449b26f4c";
    sha256 = "sha256-GbJm+WQZYjJ+dv5O9fFC5pfhE68ZA11PyFY/E4mRiks=";
  };

  buildInputs = [
    pkgs.cmake
  ];

  configurePhase = ''
    cmake . -DCMAKE_POLICY_VERSION_MINIMUM=3.5
  '';

  buildPhase = ''
    make
  '';

  installPhase = ''
    mkdir -p $out/bin
    mv bin/seqToProfile  bin/simuReads $out/bin
  '';
}
