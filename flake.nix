{
  description = "Description for the project";

  inputs = {
    flake-parts.url = "github:hercules-ci/flake-parts";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs =
    inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
        "x86_64-darwin"
      ];
      perSystem =
        { pkgs, ... }:
        {
          devShells.default =
            let
              jupyterBase =
                ps: with ps; [
                  jupyter
                  ipykernel

                  #For magma nvim
                  pynvim
                  jupyter-client
                  ueberzug
                  pillow
                  cairosvg
                  pnglatex
                  plotly
                  pyperclip

                  jupyterlab-lsp
                  python-lsp-server
                  jupyterlab-vim
                ];

              physicsPackages =
                ps:
                (
                  (jupyterBase ps)
                  ++ (with ps; [
                    numpy
                    matplotlib
                    scipy
                    pandas
                    astropy
                  ])
                );

            in
            pkgs.mkShell {
              packages = with pkgs; [
                root
                cmake
                pkg-config
                (python3.withPackages physicsPackages)
              ];
              shellHook = ''
                printf 'CompileFlags:\n  Add: [-std=c++17, %s]\n' \
                  "$(root-config --cflags | sed 's/ /, /g')" > .clangd

                mkdir -p /tmp/.jupyter
                export JUPYTER_CONFIG_DIR=/tmp/.jupyter
                python -m ipykernel install --user \
                --name physics-jupyter \
                --display-name "Nix Jupyter Environment"
              '';
            };
          packages.scripts = pkgs.stdenv.mkDerivation {
            src = ./.;
            pname = "H4l-scripts";
            version = "0.1";
            buildInputs = [ pkgs.root ];
            nativeBuildInputs = with pkgs; [
              cmake
              pkg-config
            ];
          };
        };
    };
}
