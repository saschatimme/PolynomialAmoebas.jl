using Documenter, PolynomialAmoebas

makedocs(
    sitename = "PolynomialAmoebas.jl",
    pages = [
        "Introduction" => "index.md",
        "Amoeba" => "amoeba.md",
        "Spine" => "spine.md",
        "Contour" => "contour.md",
        "Coamoeba" => "coamoeba.md",
        "Imaginary Projection" => "imaginary.md",
        "Reference" => "reference.md"
    ]
)


deploydocs(
    repo = "github.com/saschatimme/PolynomialAmoebas.jl.git",
)
