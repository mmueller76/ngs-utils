package org.nxtgenutils.util;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.net.URL;

/**
 * This file is part of NxtGenUtils.
 *
 * NxtGenUtils is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NxtGenUtils is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NxtGenUtils.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 21-Oct-2011
 * Time: 13:52:21
 */

/**
 * This class allows the determination of file types.
 *
 * @author Lennart Martens
 */
public class FileTypeDeterminator {

    /**
     * This method will try to determine whether this stream is connected to a zip file,
     * a gzipped file, a possible text file or a possible binary file, in that order.
     * It will return an appropriate Reader for the identified filetype, or 'null' if the type
     * was not recognizable (ie: we could not read the four bytes constituting the magic number).
     *
     * @param aStream   InputStream to determine the nature of
     * @return Reader   with a Reader into the converted stream.
     */
    public static Reader determineFileType(InputStream aStream) throws IOException {
        Reader result = null;
        try {
            // Reading the first four bytes.
            byte[] firstFour = new byte[4];
            PushbackInputStream pbi = new PushbackInputStream(aStream, 4);
            int count = pbi.read(firstFour, 0, 4);
            // If we couldn't even read 4 bytes, it cannot be good. In that case we'll therefore return 'null'.
            if(count == 4) {
                //logger.debug("Read magic numbers (first four bytes: " + firstFour[0] + " " + firstFour[1] + " " + firstFour[2] + " " + firstFour[3] + " " + ").");
                // Now unread our bytes.
                pbi.unread(firstFour);
                // Okay, let's check these magic numbers, shall we?
                if(firstFour[0] == (byte)0x1F && firstFour[1] == (byte)0x8b) {
                    // GZIP!
                    //logger.debug("Detected GZIP format.");
                    result = new InputStreamReader(new GZIPInputStream(pbi));
                } else if(firstFour[0] == (byte)0x50 && firstFour[1] == (byte)0x4b && firstFour[2] == (byte)0x03 && firstFour[3] == (byte)0x04) {
                    // (pk)ZIP!
                    //logger.debug("Detected ZIP format.");
                    ZipInputStream zis = new ZipInputStream(pbi);
                    ZipEntry ze = zis.getNextEntry();
                    // Extra check: ze cannot be 'null'!
                    if(ze != null) {
                        result = new InputStreamReader(zis);
                    }
                }
                // If we are here and result is still 'null', we weren't able to identify it as either GZIP or ZIP.
                // So create a regular reader and go from here.
                if(result == null) {
                    //logger.debug("Defaulted to standard Reader.");
                    result = new InputStreamReader(pbi);
                }
            }
        } catch(IOException ioe) {
            throw ioe;
            //logger.error("IOException while attempting to determine filetype: " + ioe.getMessage());
        }

        return result;
    }

    public static Reader determineFileType(URL aURL) throws IOException {

            InputStream is = null;

            //open stream according to protocol
            if (aURL.getProtocol().equals("file")) is = aURL.openStream();
            else if (aURL.getProtocol().equals("ftp")) is = aURL.openConnection().getInputStream();
            else if (aURL.getProtocol().equals("http")) is = aURL.openConnection().getInputStream();
            else throw new IOException("Unsupported protocol.");

            //determine file format
            return FileTypeDeterminator.determineFileType(is);

    }

}
