package nhs.genetics.cardiff.framework;

import java.util.IllegalFormatException;

/**
 * Calculate Hamming (edit) Distance
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2016-10-18
 */
public class Hamming
{
    public static int getHammingDistance(String s, String t) throws IllegalArgumentException
    {
        int n = 0;

        if (s.length() != t.length())
        {
            throw new IllegalArgumentException("Strings must be the same length");
        }

        for (int i = 0; i < s.length(); i++)
        {
            if (s.charAt(i) != t.charAt(i)) {
                n++;
            }
        }

        return n;
    }
}