{
  inputs.nixpkgs.url = "nixpkgs/nixos-20.03";

  outputs = { self, nixpkgs }: {

    devShell.x86_64-linux =
      with import nixpkgs { system = "x86_64-linux"; };
      stdenv.mkDerivation {
        name = "dhdm";
        buildInputs = [ mesa_glu glew glfw libpng glm fmt ];
      };

  };

}
