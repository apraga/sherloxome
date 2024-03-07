{

    inputs = {
        nixpkgs.url = "nixpkgs/nixos-23.11";
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
            };
        };
}
