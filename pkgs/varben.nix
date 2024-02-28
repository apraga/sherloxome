{ pkgs, stdenv, fetchFromGitHub, coreutils, boost, zlib, lib, makeWrapper }:

let
  # Bcftools needs perl
  runtimeInputs  = with pkgs; [
    bwa
    coreutils
    bcftools
    bedtools
    samtools
  ];
  my-python-packages = p: with p; [
    pysam
    numpy
  ];
 my-python = pkgs.python3.withPackages my-python-packages;
in
stdenv.mkDerivation {
  name = "varben";
  src = fetchFromGitHub {
    owner = "nccl-jmli";
    repo = "VarBen";
    rev = "0f66e35dc85b80938df2beafa5919330c9356953";
    sha256 ="sha256-hTlg5w4kyfsEgAINlrRxW3ZLbZ4yzm1VUNt56HmNazU=";
  };
  patches = [
    ./varben-python3.patch
];
  buildInputs = with pkgs; [
    my-python
  ];
  nativeBuildInputs = [ pkgs.makeWrapper ];

  propagatedBuildInputs = runtimeInputs ;
  installPhase = ''
  mkdir $out
  cp -r bin src varben $out/
  '';
  postFixup = ''
    makeWrapper ${my-python}/bin/python $out/bin/muteditor \
        --set PATH ${lib.makeBinPath runtimeInputs }  \
        --add-flags "$out/bin/muteditor.py"
    '';
}
