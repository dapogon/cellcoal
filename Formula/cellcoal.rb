class Cellcoal < Formula
  desc "Coalescent simulator of single-cell NGS data (somatic mutations, ADO, doublets, reads)"
  homepage "https://github.com/dapogon/cellcoal"
  url "https://github.com/dapogon/cellcoal/archive/refs/tags/v2.2.6.tar.gz"
  sha256 "450341260d0746364134af68590621af695782cd34251a2dcc9d55f6778ec0fd"
  license "GPL-2.0-or-later"

  depends_on "libomp"

  def install
    system "make"
    system "make", "install", "PREFIX=#{prefix}"
  end

  test do
    system bin/"cellcoal", "-n1", "-s10", "-l200", "-u1e-6", "-e10000", "-b1", "-#1"
    assert_predicate testpath/"Results", :directory?
  end
end
