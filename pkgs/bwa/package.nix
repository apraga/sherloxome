{ lib, stdenv, fetchFromGitHub, zlib }:

stdenv.mkDerivation {
  pname = "bwa";
  version = "0.7.17-r1198-dirty";

  src = fetchFromGitHub {
    owner = "lh3";
    repo = "bwa";
    rev = "139f68fc4c3747813783a488aef2adc86626b01b";
    hash = "sha256-8u35lTK6gBKeapYoIkG9MuJ/pyy/HFA2OiPn+Ml2C6c=";
  };

  buildInputs = [ zlib ];

  preConfigure = ''
    sed -i '/^CC/d' Makefile
  '';

  installPhase = ''
    runHook preInstall
    install -vD -t $out/bin bwa
    install -vD -t $out/lib libbwa.a
    install -vD -t $out/include bntseq.h
    install -vD -t $out/include bwa.h
    install -vD -t $out/include bwamem.h
    install -vD -t $out/include bwt.h
    runHook postInstall
  '';

  meta = {
    description = "Burrows-Wheeler Aligner for short-read alignment";
    license = lib.licenses.gpl3Plus;
    platforms = lib.platforms.linux;
  };
}
