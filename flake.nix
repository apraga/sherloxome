{

    inputs = {
        nixpkgs.url = "nixpkgs/nixos-25.11";
    };

    description = "Exome validator";

    outputs = { self, nixpkgs }:
        let
            system = "x86_64-linux";
            pkgs = import nixpkgs { inherit system; };
        in {
            packages.${system} = {
                varben = pkgs.callPackage pkgs/varben.nix {};
                simuscop = pkgs.callPackage pkgs/simuscop.nix {};
                bwa = pkgs.bwa;
                hap-py = pkgs.hap-py;
                samtools = pkgs.samtools;
            };
        };
}
