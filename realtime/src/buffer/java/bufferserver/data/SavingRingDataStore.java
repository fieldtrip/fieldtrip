package nl.fcdonders.fieldtrip.bufferserver.data;

import nl.fcdonders.fieldtrip.bufferserver.exceptions.DataException;
import nl.fcdonders.fieldtrip.bufferserver.network.NetworkProtocol;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class SavingRingDataStore extends RingDataStore {
    private static final int dataBufSize = 1024 * 10 * 10; // ~10s data
    private static final int eventBufSize = 100 * 10 * 10;  // ~10s event file-write buffer
    private BufferedOutputStream eventWriter;
    private BufferedOutputStream dataWriter;
    private BufferedOutputStream headerWriter;
    private String savePath;
    private int numReset = 0;
    private ByteBuffer writeBuf;

    /**
     * Constructor
     *
     * @param nBuffer Capacity of event buffer
     */
    public SavingRingDataStore(final int nBuffer) {
        super(nBuffer);
        initFiles();
    }


    public SavingRingDataStore(final int nBuffer, final String path) {
        super(nBuffer);
        initFiles(path);
    }

    public SavingRingDataStore(final int nBuffer, final File theDir) {
        super(nBuffer);
        initFiles(theDir);
    }


    /**
     * Constructor
     *
     * @param nSamples Capacity of the sample ringbuffer.
     * @param nEvents  Capacity of the event ringbuffer.
     */
    public SavingRingDataStore(final int nSamples, final int nEvents) {
        super(nSamples, nEvents);
        initFiles();
    }

    public SavingRingDataStore(final int nSamples, final int nEvents, String path) {
        super(nSamples, nEvents);
        initFiles(path);
    }

    public SavingRingDataStore(final int nSamples, final int nEvents, File theDir) {
        super(nSamples, nEvents);
        initFiles(theDir);
    }

    void initFiles() {
        initFiles(".");
    }

    void initFiles(String path) {
        if (path == null) path = ".";
        // add the reset number prefix to the path
        String fullpath = path + File.separator + String.format("%03d", numReset);
        File theDir = new File(fullpath);
    }

    void initFiles(File file) {
        try {
            // if the directory does not exist, create it
            if (!file.exists()) {
                try {
                    file.mkdirs();
                } catch (SecurityException se) {
                    //handle it
                }
            }
            // record the save path used
            savePath = file.getPath();
            dataWriter = new BufferedOutputStream(new FileOutputStream(savePath + File.separator + "samples"),
                    dataBufSize);
            eventWriter = new BufferedOutputStream(new FileOutputStream(savePath + File.separator + "events"),
                    eventBufSize);
            headerWriter = new BufferedOutputStream(new FileOutputStream(savePath + File.separator + "header"));
            // Write everything in BIG_ENDIAN
            writeBuf = ByteBuffer.allocate(eventBufSize);
            writeBuf.order(ByteOrder.nativeOrder());
        } catch (IOException x) {
            System.err.println(x);
        }
    }

    void resetFiles() {
        cleanup();
        initFiles(savePath);
    }

    /**
     * Removes all data.
     *
     * @throws DataException
     */
    @Override
    public synchronized void flushData() throws DataException {
        super.flushData();
        resetFiles(); // resets the sample counter, so need start new save file
    }

    /**
     * Removes the header, and all data & events.
     *
     * @throws DataException
     */
    @Override
    public synchronized void flushHeader() throws DataException {
        super.flushHeader();
        resetFiles(); // resets the sample counter, so need start new save file
    }


    /**
     * Appends the data to the storage. Throws DataException if impossible.
     *
     * @param data
     * @throws DataException
     */
    @Override
    public synchronized int putData(final Data data) throws DataException {
        if (data.dataType != dataType) {
            throw new DataException("Trying to append data of wrong dataType.");
        }
        if (data.nChans != nChans) {
            throw new DataException("Trying to append data with wrong number of channels");
        }

        // Check if byte order needs to be flipped
        if (data.order != NATIVE_ORDER && nBytes != 1) {
            for (int i = 0; i < data.nSamples; i++) {
                final byte[][] sample = new byte[nChans][nBytes];

                for (int j = 0; j < nChans; j++) {
                    for (int k = 0; k < nBytes; k++) {
                        sample[j][k] = data.data[i][j][nBytes - k - 1];
                    }
                }

                dataBuffer.add(sample);
                try {
                    dataWriterwrite(sample);
                } catch (IOException e) {
                    System.err.println("IOException writing data");
                }
            }
        } else {
            for (int i = 0; i < data.nSamples; i++) {
                dataBuffer.add(data.data[i]);
                try {
                    dataWriterwrite(data.data[i]);
                } catch (IOException e) {
                    System.err.println("IOException writing data");
                }
            }
        }
        checkListeners();
        return dataBuffer.sampleCount();
    }

    /**
     * Appends the events to the storage. Throws DataException if impossible.
     *
     * @param events
     * @throws DataException
     */
    @Override
    public synchronized int putEvents(final Event[] events) throws DataException {
        for (final Event event : events) {
            if (event.order != NATIVE_ORDER) {
                final int typeNBytes = NetworkProtocol.dataTypeSize(event.typeType);

                final byte[][] type = event.type.clone();
                if (typeNBytes > 1) {
                    for (int i = 0; i < event.typeSize; i++) {
                        for (int j = 0; j < typeNBytes; j++) {
                            type[i * typeNBytes + j] = event.type[i * typeNBytes + typeNBytes - j - 1];
                        }
                    }
                }

                final int valueNBytes = NetworkProtocol.dataTypeSize(event.valueType);

                final byte[][] value = event.value.clone();
                if (valueNBytes > 1) {
                    for (int i = 0; i < event.valueSize; i++) {
                        for (int j = 0; j < valueNBytes; j++) {
                            value[i * valueNBytes + j] = event.value[i * valueNBytes + valueNBytes - j - 1];
                        }
                    }
                }

                Event evt = new Event(event, type, value, NATIVE_ORDER);
                eventBuffer.add(evt);
                try {
                    eventWriterwrite(evt);
                } catch (IOException e) {
                    System.err.println("IOException writing events");
                }
            } else {
                eventBuffer.add(event);
                try {
                    eventWriterwrite(event);
                } catch (IOException e) {
                    System.err.println("IOException writing events");
                }
            }
        }
        checkListeners();
        return eventBuffer.eventCount();
    }

    /**
     * Adds the header to the storage. Throws DataException if impossible.
     *
     * @param hdr
     * @throws DataException
     */
    @Override
    public synchronized void putHeader(Header header) throws DataException {

        final boolean newHeader = header == null;

        // Check if header is in BIG_ENDIAN ByteOrder.
        if (header.order != NATIVE_ORDER) {
            final Chunk[] chunks = header.chunks;

            // Check each chunk, if it is a CHUNK_RESOLUTIONS chunk, flip the
            // byte order.
            for (int i = 0; i < chunks.length; i++) {
                if (chunks[i].type == NetworkProtocol.CHUNK_RESOLUTIONS) {
                    final byte[] data = new byte[chunks[i].data.length];

                    for (int j = 0; j < header.nChans; j++) {
                        for (int k = 0; k < 8; k++) {
                            data[j * 8 + k] = chunks[i].data[j * 8 + 7 - k];
                        }
                    }

                    // Replace chunk.
                    chunks[i] = new Chunk(chunks[i].type, chunks[i].size, data);
                }
            }

            // Create new header with BIG_ENDIAN ByteOrder
            header = new Header(header, chunks, NATIVE_ORDER);
        }

        if (newHeader) {
            if (nChans != header.nChans) {
                throw new DataException("Replacing header has different number of channels");
            }
            if (dataType != header.dataType) {
                throw new DataException("Replacing header has different data type");
            }
        } else {
            nChans = header.nChans;
            dataType = header.dataType;
            nBytes = NetworkProtocol.dataTypeSize(dataType);
        }

        this.header = header;
        dataBuffer = new DataRingBuffer(dataBufferSize, nChans, nBytes);
        try {
            headerWriterwrite(header);
        } catch (IOException e) {
            System.err.println("IOException writing header");
        }
    }


    // Methods to write to save files
    void dataWriterwrite(byte[][] sample) throws IOException {
        if (dataWriter == null) return;
        //System.err.println("Writing samples:");
        int n = 0;
        for (int i = 0; i < sample.length; i++) {
            dataWriter.write(sample[i]);
            n += sample[i].length;
        }
        //System.err.println("wrote"+n+"bytes");
        //dataWriter.flush();
    }

    void eventWriterwrite(Event event) throws IOException {
        if (eventWriter == null) return;
        //System.err.println("Writing event:");
        writeBuf.clear();
        event.serialize(writeBuf);
        eventWriter.write(writeBuf.array(), writeBuf.arrayOffset(), writeBuf.position());
        //System.err.println("Event len"+writeBuf.position()+"char");
        //eventWriter.flush();
    }

    void headerWriterwrite(Header header) throws IOException {
        if (headerWriter == null) return;
        //System.err.println("Writing header:");
        writeBuf.clear();
        header.serialize(writeBuf);
        headerWriter.write(writeBuf.array(), writeBuf.arrayOffset(), writeBuf.position());
        //System.err.println(writeBuf.position()+"char");
        headerWriter.flush(); // force flush header... we always want this correct
    }

    public void cleanup() {
        if (headerWriter != null) try {
            headerWriter.close();
        } catch (IOException e) {
            System.err.println("IOException closing header");
        }

        if (eventWriter != null) try {
            eventWriter.close();
        } catch (IOException e) {
            System.err.println("IOException closing events");
        }
        if (dataWriter != null) try {
            dataWriter.close();
        } catch (IOException e) {
            System.err.println("IOException closing samples");
        }
    }
}
