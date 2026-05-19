{ lib, stdenv, fetchFromGitHub, zlib }:

stdenv.mkDerivation {
  pname = "bwa";
  version = "0.7.18";

  src = fetchFromGitHub {
    owner = "lh3";
    repo = "bwa";
    rev = "v0.7.18";
    hash = "sha256-ITvugdgUUfncDcJjEcBaO8ux2fZ4YPEdg3/i/iePw+0=";
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
