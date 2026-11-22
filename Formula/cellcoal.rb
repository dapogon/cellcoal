class Cellcoal < Formula
  desc "Coalescent simulator of single-cell NGS data (somatic mutations, ADO, doublets, reads)"
  homepage "https://github.com/dapogon/cellcoal"
  url "https://github.com/dapogon/cellcoal/archive/refs/tags/v2.6.0.tar.gz"
  sha256 "3a631ac8f6b3fe2d383695c2ac3e083ee347ad71b73d862575dfb37968e8e18c"
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
