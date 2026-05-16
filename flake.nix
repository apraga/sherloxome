{

    inputs = {
        nixpkgs.url = "nixpkgs/nixos-25.11";
    };

    description = "Exome validator";

    outputs = { self, nixpkgs }:
        let
            system = "x86_64-linux";
            pkgs = import nixpkgs { inherit system; };
            local-rtg-tools = pkgs.callPackage pkgs/rtg-tools/package.nix {};
        in {
            packages.${system} = {
                varben = pkgs.callPackage pkgs/varben/package.nix {};
                simuscop = pkgs.callPackage pkgs/simuscop/package.nix {};
                bwa = pkgs.bwa;
                # Wait for upstream to add it
                # rtg-tools = pkgs.rtg-tools;
                # hap-py = pkgs.hap-py;
                rtg-tools = local-rtg-tools;
                hap-py = pkgs.callPackage pkgs/hap-py/package.nix { rtg-tools = local-rtg-tools; };
                samtools = pkgs.samtools;
            };
        };
}
