{
  description = "Development environment for TDA362 Computer Graphics Labs";

  inputs = {
    # Use unstable for most packages
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    # Pin nixos-24.05 for embree2
    nixpkgs-24-05.url = "github:NixOS/nixpkgs/nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, nixpkgs-24-05, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        # Import nixpkgs-24.05 for embree2
        pkgs-24-05 = nixpkgs-24-05.legacyPackages.${system};

        # Import unstable with overlay to add embree2 from 24.05
        pkgs = import nixpkgs {
          inherit system;
          overlays = [
            (final: prev: {
              embree2 = pkgs-24-05.embree2;
            })
          ];
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            # Build tools
            cmake
            gnumake
            gcc

            # Graphics libraries
            SDL2
            glew
            glm
            mesa
            libGL
            libGLU

            # X11 libraries (needed by SDL2 headers)
            xorg.libX11
            xorg.libXext
            xorg.libXcursor
            xorg.libXi
            xorg.libXrandr
            xorg.libXxf86vm

            # Pathtracer dependency (embree 2.x from nixos-24.05)
            embree2

            # Additional utilities
            pkg-config
          ];

          # Set environment variables for CMake to find libraries
          CMAKE_PREFIX_PATH = "${pkgs.SDL2}:${pkgs.glew}:${pkgs.glm}:${pkgs.embree2}";

          # Ensure OpenGL libraries are found
          LD_LIBRARY_PATH = "${pkgs.lib.makeLibraryPath [
            pkgs.libGL
            pkgs.libGLU
            pkgs.SDL2
            pkgs.glew
            pkgs.embree2
          ]}";
        };
      }
    );
}

