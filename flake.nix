{

    inputs = {
        nixpkgs.url = "github:NixOS/nixpkgs/68799d2b3b33ca5a2406d99ec5d67900a1ea658f";
        # nixpkgs.url = "nixpkgs/nixos-25.11";
    };

    description = "Exome validator";

    outputs = { self, nixpkgs }:
        let
            system = "x86_64-linux";
            pkgs = import nixpkgs { inherit system; };
            deps = {
                varben = pkgs.callPackage pkgs/varben/package.nix {};
                simuscop = pkgs.callPackage pkgs/simuscop/package.nix {};
                bwa = pkgs.bwa;
                rtg-tools = pkgs.rtg-tools;
                hap-py = pkgs.hap-py;
                samtools = pkgs.samtools;
            };
        in {
            packages.${system} = deps // {
                default = pkgs.buildEnv {
                    name = "sherloxome";
                    paths = builtins.attrValues deps;
                };
            };
        };
}
