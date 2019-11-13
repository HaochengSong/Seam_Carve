import edu.princeton.cs.algs4.Picture;
import java.awt.Color;
public class SeamCarver {
	private Picture picture;
	private double[][] energy;
	public SeamCarver(Picture picture) {
		if (picture == null)
			throw new IllegalArgumentException("Illeagal!!");
		this.picture = new Picture(picture);
		energy = new double[height()][width()];
		for (int i = 0; i < height(); i++)
			for (int j = 0; j < width(); j++)
				energy[i][j] = computeEnergy(j, i);
	}
	private double calX(int x, int y) {
		Color c1 = picture.get(x - 1, y);
		Color c2 = picture.get(x + 1, y);
		int r = c1.getRed() - c2.getRed();
		int g = c1.getGreen() - c2.getGreen();
		int b = c1.getBlue() - c2.getBlue();
		return r * r + g * g + b * b;
	}
	private double calY(int x, int y) {
		Color c1 = picture.get(x, y - 1);
		Color c2 = picture.get(x, y + 1);
		int r = c1.getRed() - c2.getRed();
		int g = c1.getGreen() - c2.getGreen();
		int b = c1.getBlue() - c2.getBlue();
		return r * r + g * g + b * b;
	}
	private double computeEnergy(int x, int y) {
        if (x == 0 || x == this.width() - 1 || y == 0 || y == this.height() - 1)
            return 1000.0;
        else
        	return Math.sqrt(calX(x ,y) + calY(x, y));
	}
	public Picture picture() {
		return new Picture(this.picture);
	}
	public int width() {
		return picture.width();
	}
	public int height() {
		return picture.height();
	}
	public double energy(int x, int y) {
		if (x < 0 || x > this.width() - 1 || y < 0 || y > this.height() - 1) 
            throw new IllegalArgumentException("Illeagal!!");
		return energy[y][x];
	}
	public int[] findHorizontalSeam() {
		double[][] distTo = new double[height()][width()];
		int[][] edgeTo = new int[height()][width()];
		for (int i = 0; i < height(); i++)
			for (int j = 1; j < width(); j++)
				distTo[i][j] = Double.MAX_VALUE;
		for (int column = 0; column < width() - 1; column++)
			for (int row = 1; row < height() - 1; row++) {
				for (int i = -1; i < 2; i++) {
					if (distTo[row][column] + energy[row + i][column + 1] < distTo[row + i][column + 1]) {
						distTo[row + i][column + 1] = distTo[row][column] + energy[row + i][column + 1];
						edgeTo[row + i][column + 1] = row;
					}
				}
			}
		double minPath = Double.MAX_VALUE;
		int minRow = 0;
		for (int i = 0; i < height(); i++)
			if (distTo[i][width() - 1] < minPath) {
				minPath = distTo[i][width() - 1];
				minRow = i;
			}
		int[] seam = new int[width()];
		seam[width() - 1] = minRow;
		for (int i = width() - 2; i >= 0; i--) {
			seam[i] = edgeTo[minRow][i + 1];
			minRow = edgeTo[minRow][i + 1];
		}
		return seam;
	}
	public int[] findVerticalSeam() {
		double[][] distTo = new double[height()][width()];
		int[][] edgeTo = new int[height()][width()];
		for (int i = 1; i < height(); i++)
			for (int j = 0; j < width(); j++)
				distTo[i][j] = Double.MAX_VALUE;
		for (int row = 0; row < height() - 1; row++)
			for (int column = 1; column < width() - 1; column++) {
				for (int i = -1; i < 2; i++) {
					if (distTo[row][column] + energy[row + 1][column + i] < distTo[row + 1][column + i]) {
						distTo[row + 1][column + i] = distTo[row][column] + energy[row + 1][column + i];
						edgeTo[row + 1][column + i] = column;
					}			
				}
			}
		double minPath = Double.MAX_VALUE;
		int minColumn = 0;
		for (int i = 0; i < width(); i++)
			if (distTo[height() - 1][i] < minPath) {
				minPath = distTo[height() - 1][i];
				minColumn = i;
			}
		int[] seam = new int[height()];
		seam[height() - 1] = minColumn;
		for (int i = height() - 2; i >= 0; i--) {
			seam[i] = edgeTo[i + 1][minColumn];
			minColumn = edgeTo[i + 1][minColumn];
		}
		return seam;
	}
	public void removeHorizontalSeam(int[] seam) {
        if (height() <= 1) throw new IllegalArgumentException("Illeagal!!");
        if (seam == null) throw new IllegalArgumentException("Illeagal!!");
        if (seam.length != width()) throw new IllegalArgumentException("Illeagal!!");
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] > height() - 1)
                throw new IllegalArgumentException();
            if (i < width() - 1 && Math.pow(seam[i] - seam[i + 1], 2) > 1)
                throw new IllegalArgumentException();
        }
		
		Picture newPicture = new Picture(width(), height() - 1);
		for (int column = 0; column < width(); column++) {
			for (int row = 0; row < height() - 1; row++) {
				if (row < seam[column]) {
					newPicture.setRGB(column, row, picture.getRGB(column, row));
				}
				else {
					newPicture.setRGB(column, row, picture.getRGB(column, row + 1));
				}
			}
		}
		this.picture = newPicture;
		double[][] newEnergy = new double[height()][width()];
		for (int column = 0, row = 0; column < width(); column++) {
			row = 0;
			int index = seam[column];
			if (index == 0) {
				newEnergy[row][column] = computeEnergy(column, row);
				row++;
				for (; row < height(); row++)
					newEnergy[row][column] = energy[row + 1][column];
			}
			else if (index == height()) {
				for (; row < index - 1; row++)
					newEnergy[row][column] = energy[row][column];
				newEnergy[row][column] = computeEnergy(column, row);
			}
			else {
				for (; row < index - 1; row++)
					newEnergy[row][column] = energy[row][column];
				newEnergy[row][column] = computeEnergy(column, row);
				row++;
				newEnergy[row][column] = computeEnergy(column, row);
				row++;
				for (; row < height(); row++)
					newEnergy[row][column] = energy[row + 1][column];
			}				
		}
		this.energy = newEnergy;
	}
	public void removeVerticalSeam(int[] seam) {
        if (width() <= 1) throw new IllegalArgumentException("Illeagal!!");
        if (seam == null) throw new IllegalArgumentException("Illeagal!!");
        if (seam.length != height()) throw new IllegalArgumentException("Illeagal!!");
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] > width() - 1)
                throw new IllegalArgumentException("Illeagal!!");
            if (i < height() - 1 && Math.pow(seam[i] - seam[i + 1], 2) > 1)
                throw new IllegalArgumentException("Illeagal!!");
        }
		Picture newPicture = new Picture(width() - 1, height());
		for (int row = 0; row < height(); row++) {
			for (int column = 0; column < width() - 1; column++) {
				if (column < seam[row])
					newPicture.setRGB(column, row, picture.getRGB(column, row));
				else
					newPicture.setRGB(column, row, picture.getRGB(column + 1, row));
			}
		}
		this.picture = newPicture;
		double[][] newEnergy = new double[height()][width()];
		for (int row = 0, column = 0; row < height(); row++) {
			column = 0;
			int index = seam[row];
			if (index == 0) {
				newEnergy[row][column] = computeEnergy(column, row);
				column++;
				for (; column < width(); column++)
					newEnergy[row][column] = energy[row][column + 1];
			}
			else if (index == width()) {
				for (; column < index - 1; column++)
					newEnergy[row][column] = energy[row][column];
				newEnergy[row][column] = computeEnergy(column, row);
			}
			else {
				for (; column < index - 1; column++)
					newEnergy[row][column] = energy[row][column];
				newEnergy[row][column] = computeEnergy(column, row);
				column++;
				newEnergy[row][column] = computeEnergy(column, row);
				column++;
				for (; column < width(); column++)
					newEnergy[row][column] = energy[row][column + 1];
			}
		}
		this.energy = newEnergy;
	}
    /*public static void main(String[] args) {
        if (args.length != 3) {
            StdOut.println("Usage:\njava ResizeDemo [image filename] [num cols to remove] [num rows to remove]");
            return;
        }

        Picture inputImg = new Picture(args[0]);
        int removeColumns = Integer.parseInt(args[1]);
        int removeRows = Integer.parseInt(args[2]); 

        StdOut.printf("image is %d columns by %d rows\n", inputImg.width(), inputImg.height());
        SeamCarver sc = new SeamCarver(inputImg);

        Stopwatch sw = new Stopwatch();

        for (int i = 0; i < removeRows; i++) {
            int[] horizontalSeam = sc.findHorizontalSeam();
            sc.removeHorizontalSeam(horizontalSeam);
        }

        for (int i = 0; i < removeColumns; i++) {
            int[] verticalSeam = sc.findVerticalSeam();
            sc.removeVerticalSeam(verticalSeam);
        }
        Picture outputImg = sc.picture();

        StdOut.printf("new image size is %d columns by %d rows\n", sc.width(), sc.height());

        StdOut.println("Resizing time: " + sw.elapsedTime() + " seconds.");
        inputImg.show();
        outputImg.show();
        StdOut.println("width: " + sc.width() + " height: " + sc.height());
    }
    public static void main(String[] args) {
        Picture picture = new Picture(args[0]);
        StdOut.printf("image is %d pixels wide by %d pixels high.\n", picture.width(), picture.height());
        
        SeamCarver sc = new SeamCarver(picture);
        
        StdOut.printf("Printing energy calculated for each pixel.\n");        

        for (int row = 0; row < sc.height(); row++) {
            for (int col = 0; col < sc.width(); col++)
                StdOut.printf("%9.0f ", sc.energy(col, row));
            StdOut.println();
        }
        int[] a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        int[] b = {0, 1, 2, 3, 4, 3, 4, 5, 5};
        sc.removeHorizontalSeam(a);
        sc.removeVerticalSeam(b);
        StdOut.printf("Printing energy calculated for each pixel after removing.\n");        

        for (int row = 0; row < sc.height(); row++) {
            for (int col = 0; col < sc.width(); col++)
                StdOut.printf("%9.0f ", sc.energy(col, row));
            StdOut.println();
        }
    }*/
}
