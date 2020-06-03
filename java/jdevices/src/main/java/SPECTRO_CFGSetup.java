
import javax.swing.JLabel;
import javax.swing.table.TableModel;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author taliercio
 */
public class SPECTRO_CFGSetup extends DeviceSetup
{
	private static final long serialVersionUID = 1L;
    /**
     * Creates new form SPECTRO_CFGSetup
     */
    public SPECTRO_CFGSetup() {
        initComponents();        
    }
    public static void main(String[] args) {
        SPECTRO_CFGSetup spc = new SPECTRO_CFGSetup();
        spc.pack();
        spc.setVisible(true);
    }
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel5 = new javax.swing.JPanel();
        jPanel1 = new javax.swing.JPanel();
        deviceField4 = new DeviceField();
        jPanel2 = new javax.swing.JPanel();
        deviceField2 = new DeviceField();
        deviceField3 = new DeviceField();
        jScrollPane2 = new javax.swing.JScrollPane();
        jTable1 = new javax.swing.JTable();
        jPanel3 = new javax.swing.JPanel();
        jPanel4 = new javax.swing.JPanel();
        deviceOk2 = new DeviceOk();
        deviceApply1 = new DeviceApply();
        deviceReset1 = new DeviceReset();
        deviceCancel1 = new DeviceCancel();

        setDeviceProvider("scpsl.nbtf");
        setDeviceTitle("SPECTROSCOPY Diagnostic Configuration Device");
        setDeviceType("SPECTRO_CFG");
        setHeight(1000);
        setWidth(1300);

        jPanel5.setLayout(new java.awt.GridLayout(2, 0));

        deviceField4.setIdentifier("");
        deviceField4.setLabelString("Comment:");
        deviceField4.setNumCols(80);
        deviceField4.setOffsetNid(1);
        deviceField4.setTextOnly(true);
        jPanel1.add(deviceField4);

        jPanel5.add(jPanel1);

        deviceField2.setIdentifier("");
        deviceField2.setLabelString("Configuration File:");
        deviceField2.setNumCols(55);
        deviceField2.setOffsetNid(2);
        deviceField2.setTextOnly(true);
        jPanel2.add(deviceField2);

        deviceField3.setIdentifier("");
        deviceField3.setLabelString("Date:");
        deviceField3.setNumCols(20);
        deviceField3.setOffsetNid(3);
        deviceField3.setTextOnly(true);
        jPanel2.add(deviceField3);

        jPanel5.add(jPanel2);

        getContentPane().add(jPanel5, java.awt.BorderLayout.NORTH);

        jTable1.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null},
                {null, null, null, null, null, null, null, null, null, null, null, null, null, null, null}
            },
            new String [] {
                "#", "Label Los", "Telescope", "c-patch", "b-patch", "or-pannel[rc]", "a-patch", "Diaphragm", "Diameter", "P0", "P1", "RoI", "Spec", "Note", "Calib."
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Integer.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.Boolean.class
            };
            boolean[] canEdit = new boolean [] {
                false, false, false, false, false, false, false, false, false, false, false, false, false, true, true
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        jTable1.getTableHeader().setReorderingAllowed(false);
        jTable1.addComponentListener(new java.awt.event.ComponentAdapter() {
            public void componentShown(java.awt.event.ComponentEvent evt) {
                jTable1ComponentShown(evt);
            }
        });
        jTable1.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
            public void propertyChange(java.beans.PropertyChangeEvent evt) {
                jTable1PropertyChange(evt);
            }
        });
        jScrollPane2.setViewportView(jTable1);
        javax.swing.table.DefaultTableCellRenderer rightRenderer = new javax.swing.table.DefaultTableCellRenderer();
        rightRenderer.setHorizontalAlignment(JLabel.CENTER);
        for(int i = 0; i < jTable1.getColumnCount(); i++)
        jTable1.getColumnModel().getColumn(i).setCellRenderer(rightRenderer);

        getContentPane().add(jScrollPane2, java.awt.BorderLayout.CENTER);

        jPanel3.setLayout(new java.awt.GridLayout(1, 0));

        deviceOk2.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mousePressed(java.awt.event.MouseEvent evt) {
                deviceOk2MousePressed(evt);
            }
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                deviceOk2MouseReleased(evt);
            }
        });
        jPanel4.add(deviceOk2);

        deviceApply1.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                deviceApply1MouseReleased(evt);
            }
        });
        jPanel4.add(deviceApply1);

        deviceReset1.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                deviceReset1MouseReleased(evt);
            }
        });
        jPanel4.add(deviceReset1);
        jPanel4.add(deviceCancel1);

        jPanel3.add(jPanel4);

        getContentPane().add(jPanel3, java.awt.BorderLayout.SOUTH);
    }// </editor-fold>//GEN-END:initComponents

    int SPECTRO_CFG_LOS_001_NOTE    = 32;
    int SPECTRO_CFG_LOS_OFFSET      = 16;
    int MAX_LOS                     = 180;

    private void resetTable(TableModel tbm)
    {
        String str_data;
        int k, c, r;
        Boolean b = new Boolean(Boolean.FALSE);
        
        for( k = 0, c = tbm.getColumnCount()-2; c<tbm.getColumnCount() ;c++ ,k++ ) {
            for(r=0;r<tbm.getRowCount();r++)
            {
                int nidData = baseNid+SPECTRO_CFG_LOS_001_NOTE + SPECTRO_CFG_LOS_OFFSET * r + k;
                try {    
                     str_data = subtree.getString(subtree.getDataExpr(nidData));
                }catch(Exception e) {continue;}
                if ( tbm.getColumnClass(c).isInstance(b) )
                {
                    if(str_data.equals("true"))
                        tbm.setValueAt(Boolean.TRUE, r, c);
                    else
                        tbm.setValueAt(Boolean.FALSE, r, c);
                }
                else
                    tbm.setValueAt(str_data, r, c);
            }
        }    
    }
 
    private void saveTable(TableModel tbm)
    {
        String new_value, old_value;
        int k, c, r;
        
        for(k =0,c=tbm.getColumnCount()-2;c<tbm.getColumnCount();c++,k++) {
            for(r=0;r<tbm.getRowCount();r++)
            {
               try {
                   new_value = tbm.getValueAt( r, c).toString();
               } catch(Exception e) {continue;}
               
               int nidData = baseNid+SPECTRO_CFG_LOS_001_NOTE + SPECTRO_CFG_LOS_OFFSET * r + k;
               try {    
                    old_value = subtree.getString(subtree.getDataExpr(nidData));
                    if (new_value.equals(old_value))
                        continue;
                    subtree.putDataExpr(nidData, new_value);
               }catch(Exception e) {continue;}
            }
        }    
    }
    
    private void jTable1PropertyChange(java.beans.PropertyChangeEvent evt) {//GEN-FIRST:event_jTable1PropertyChange
        // TODO add your handling code here:
        int r=0,c=0;
        int offsetNid = 5;//CONFIG_LOS offset nid
        String curr_data;
        String data[];
        Boolean los_flag[] = new Boolean[MAX_LOS];
       
        if( evt.getPropertyName().equals("tableCellEditor") )
            return;
                    
        TableModel tbm = jTable1.getModel();
        
        //Column 0 is LoS number
        for(r=0;r<tbm.getRowCount();r++)
        {
                tbm.setValueAt(r+1, r, 0);
                los_flag[r] = false;
        }
        
        //Last two colums, note ad calib flag, can be changed via setup.
        //Values are read from LOS node
        for(c=1;c<tbm.getColumnCount()-2;c++) {
            int nidData = baseNid+offsetNid + (c-1);
            try {
                curr_data = subtree.getDataExpr(nidData);
                data = subtree.getStringArray(curr_data);
            }catch(Exception e) {continue;}
            
            for(r=0;r<Math.min(tbm.getRowCount(), data.length);r++)
            {
                //Only if LoS lable is defined row values table are written
                if( c == 1 )
                {
                    if( !(data[r].trim().equals("None")) && (data[r].trim().length() != 0)  )
                        los_flag[r] = true;
                }
                if(los_flag[r])
                    tbm.setValueAt(data[r], r, c);
                else
                    tbm.setValueAt("", r, c);
            }
        }
        
        
        if (jTable1.getColumnModel().getColumnCount() > 0) {
            jTable1.getColumnModel().getColumn(0).setMinWidth(30);
            jTable1.getColumnModel().getColumn(0).setPreferredWidth(30);
            jTable1.getColumnModel().getColumn(0).setMaxWidth(30);
            jTable1.getColumnModel().getColumn(1).setMinWidth(65);
            jTable1.getColumnModel().getColumn(1).setPreferredWidth(65);
            jTable1.getColumnModel().getColumn(1).setMaxWidth(65);
            jTable1.getColumnModel().getColumn(2).setMinWidth(75);
            jTable1.getColumnModel().getColumn(2).setPreferredWidth(75);
            jTable1.getColumnModel().getColumn(2).setMaxWidth(75);
            jTable1.getColumnModel().getColumn(3).setMinWidth(60);
            jTable1.getColumnModel().getColumn(3).setPreferredWidth(60);
            jTable1.getColumnModel().getColumn(3).setMaxWidth(60);
            jTable1.getColumnModel().getColumn(4).setMinWidth(60);
            jTable1.getColumnModel().getColumn(4).setPreferredWidth(60);
            jTable1.getColumnModel().getColumn(4).setMaxWidth(60);
            jTable1.getColumnModel().getColumn(5).setMinWidth(85);
            jTable1.getColumnModel().getColumn(5).setPreferredWidth(85);
            jTable1.getColumnModel().getColumn(5).setMaxWidth(85);
            jTable1.getColumnModel().getColumn(6).setMinWidth(60);
            jTable1.getColumnModel().getColumn(6).setPreferredWidth(60);
            jTable1.getColumnModel().getColumn(6).setMaxWidth(60);
            jTable1.getColumnModel().getColumn(7).setMinWidth(75);
            jTable1.getColumnModel().getColumn(7).setPreferredWidth(75);
            jTable1.getColumnModel().getColumn(7).setMaxWidth(75);
            jTable1.getColumnModel().getColumn(8).setMinWidth(65);
            jTable1.getColumnModel().getColumn(8).setPreferredWidth(65);
            jTable1.getColumnModel().getColumn(8).setMaxWidth(65);
            jTable1.getColumnModel().getColumn(9).setMinWidth(100);
            jTable1.getColumnModel().getColumn(9).setPreferredWidth(100);
            jTable1.getColumnModel().getColumn(9).setMaxWidth(100);
            jTable1.getColumnModel().getColumn(10).setMinWidth(100);
            jTable1.getColumnModel().getColumn(10).setPreferredWidth(100);
            jTable1.getColumnModel().getColumn(10).setMaxWidth(100);
            jTable1.getColumnModel().getColumn(11).setMinWidth(33);
            jTable1.getColumnModel().getColumn(11).setPreferredWidth(33);
            jTable1.getColumnModel().getColumn(11).setMaxWidth(33);
            jTable1.getColumnModel().getColumn(12).setMinWidth(100);
            jTable1.getColumnModel().getColumn(12).setPreferredWidth(100);
            jTable1.getColumnModel().getColumn(12).setMaxWidth(100);
            jTable1.getColumnModel().getColumn(14).setMinWidth(60);
            jTable1.getColumnModel().getColumn(14).setPreferredWidth(60);
            jTable1.getColumnModel().getColumn(14).setMaxWidth(60);
        }
        
        resetTable(tbm);
    }//GEN-LAST:event_jTable1PropertyChange

    private void jTable1ComponentShown(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_jTable1ComponentShown
        // TODO add your handling code here:
       
    }//GEN-LAST:event_jTable1ComponentShown

    private void deviceApply1MouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_deviceApply1MouseReleased
        // TODO add your handling code here:
        TableModel tbm = jTable1.getModel();
        saveTable(tbm);
    }//GEN-LAST:event_deviceApply1MouseReleased

    private void deviceReset1MouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_deviceReset1MouseReleased
        // TODO add your handling code here:
        TableModel tbm = jTable1.getModel();
        resetTable(tbm);
    }//GEN-LAST:event_deviceReset1MouseReleased

    private void deviceOk2MouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_deviceOk2MouseReleased
        // TODO add your handling code here:
    }//GEN-LAST:event_deviceOk2MouseReleased

    private void deviceOk2MousePressed(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_deviceOk2MousePressed
        // TODO add your handling code here:
        TableModel tbm = jTable1.getModel();
        saveTable(tbm);
    }//GEN-LAST:event_deviceOk2MousePressed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private DeviceApply deviceApply1;
    private DeviceCancel deviceCancel1;
    private DeviceField deviceField2;
    private DeviceField deviceField3;
    private DeviceField deviceField4;
    private DeviceOk deviceOk2;
    private DeviceReset deviceReset1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JTable jTable1;
    // End of variables declaration//GEN-END:variables
}
