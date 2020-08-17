{
  inputs.nixpkgs.url = "nixpkgs/nixos-20.03";

  outputs = { self, nixpkgs }: {

    defaultPackage.x86_64-linux =
      with import nixpkgs { system = "x86_64-linux"; };
      stdenv.mkDerivation {
        name = "dhdm";
        buildInputs = [ mesa_glu glew glfw libpng glm fmt nlohmann_json opensubdiv boost ];
        src = self;
        preBuild = "cd src";
        installPhase = "mkdir -p $out/bin; cp dhdm $out/bin/";
        enableParallelBuilding = true;
      };

  };

}
