/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author mdsplus
 */
public class MARTE2_MCC118Setup extends DeviceSetup {

    /**
     * Creates new form MARTE2_MCC118Setup
     */
    public MARTE2_MCC118Setup() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        deviceButtons1 = new DeviceButtons();
        jPanel1 = new javax.swing.JPanel();
        deviceField1 = new DeviceField();
        deviceField2 = new DeviceField();
        deviceField3 = new DeviceField();
        jScrollPane2 = new javax.swing.JScrollPane();
        jPanel2 = new javax.swing.JPanel();
        jPanel3 = new javax.swing.JPanel();
        deviceField4 = new DeviceField();
        jPanel4 = new javax.swing.JPanel();
        deviceField5 = new DeviceField();
        jPanel5 = new javax.swing.JPanel();
        deviceField6 = new DeviceField();
        jPanel6 = new javax.swing.JPanel();
        deviceField7 = new DeviceField();
        jPanel7 = new javax.swing.JPanel();
        deviceField8 = new DeviceField();
        jPanel8 = new javax.swing.JPanel();
        deviceField9 = new DeviceField();
        jPanel9 = new javax.swing.JPanel();
        deviceField10 = new DeviceField();
        jPanel10 = new javax.swing.JPanel();
        deviceField11 = new DeviceField();
        jPanel11 = new javax.swing.JPanel();
        deviceField12 = new DeviceField();
        jPanel12 = new javax.swing.JPanel();
        deviceField13 = new DeviceField();
        jPanel13 = new javax.swing.JPanel();
        deviceField14 = new DeviceField();
        jPanel14 = new javax.swing.JPanel();
        deviceField15 = new DeviceField();
        jPanel15 = new javax.swing.JPanel();
        deviceField16 = new DeviceField();
        jPanel16 = new javax.swing.JPanel();
        deviceField17 = new DeviceField();
        jPanel17 = new javax.swing.JPanel();
        deviceField18 = new DeviceField();
        jPanel18 = new javax.swing.JPanel();
        deviceField19 = new DeviceField();

        setDeviceProvider("spilds.rfx.local:8100");
        setDeviceTitle("MC118 ADC");
        setDeviceType("MARTE2_MCC118");
        setHeight(300);
        setWidth(600);
        getContentPane().add(deviceButtons1, java.awt.BorderLayout.PAGE_END);

        deviceField1.setIdentifier("");
        deviceField1.setLabelString("GPIO Pin: ");
        deviceField1.setNumCols(4);
        deviceField1.setOffsetNid(10);
        jPanel1.add(deviceField1);

        deviceField2.setIdentifier("");
        deviceField2.setLabelString("GPIO Device");
        deviceField2.setOffsetNid(13);
        deviceField2.setTextOnly(true);
        jPanel1.add(deviceField2);

        deviceField3.setIdentifier("");
        deviceField3.setLabelString("Tree Writer CPU Mask");
        deviceField3.setNumCols(4);
        deviceField3.setOffsetNid(20);
        jPanel1.add(deviceField3);

        getContentPane().add(jPanel1, java.awt.BorderLayout.NORTH);

        jPanel2.setLayout(new java.awt.GridLayout(16, 1));

        deviceField4.setIdentifier("");
        deviceField4.setLabelString("Ch1 Seg. Length (0 to disable write):");
        deviceField4.setOffsetNid(26);
        jPanel3.add(deviceField4);

        jPanel2.add(jPanel3);

        deviceField5.setIdentifier("");
        deviceField5.setLabelString("Ch2 Seg. Length (0 to disable write):");
        deviceField5.setOffsetNid(34);
        jPanel4.add(deviceField5);

        jPanel2.add(jPanel4);

        deviceField6.setIdentifier("");
        deviceField6.setLabelString("Ch3 Seg. Length (0 to disable write):");
        deviceField6.setOffsetNid(42);
        jPanel5.add(deviceField6);

        jPanel2.add(jPanel5);

        deviceField7.setIdentifier("");
        deviceField7.setLabelString("Ch4 Seg. Length (0 to disable write):");
        deviceField7.setOffsetNid(50);
        jPanel6.add(deviceField7);

        jPanel2.add(jPanel6);

        deviceField8.setIdentifier("");
        deviceField8.setLabelString("Ch5 Seg. Length (0 to disable write):");
        deviceField8.setOffsetNid(58);
        jPanel7.add(deviceField8);

        jPanel2.add(jPanel7);

        deviceField9.setIdentifier("");
        deviceField9.setLabelString("Ch6 Seg. Length (0 to disable write):");
        deviceField9.setOffsetNid(66);
        jPanel8.add(deviceField9);

        jPanel2.add(jPanel8);

        deviceField10.setIdentifier("");
        deviceField10.setLabelString("Ch7 Seg. Length (0 to disable write):");
        deviceField10.setOffsetNid(74);
        jPanel9.add(deviceField10);

        jPanel2.add(jPanel9);

        deviceField11.setIdentifier("");
        deviceField11.setLabelString("Ch8 Seg. Length (0 to disable write):");
        deviceField11.setOffsetNid(82);
        jPanel10.add(deviceField11);

        jPanel2.add(jPanel10);

        deviceField12.setIdentifier("");
        deviceField12.setLabelString("Ch9 Seg. Length (0 to disable write):");
        deviceField12.setOffsetNid(90);
        jPanel11.add(deviceField12);

        jPanel2.add(jPanel11);

        deviceField13.setIdentifier("");
        deviceField13.setLabelString("Ch10 Seg. Length (0 to disable write):");
        deviceField13.setOffsetNid(98);
        jPanel12.add(deviceField13);

        jPanel2.add(jPanel12);

        deviceField14.setIdentifier("");
        deviceField14.setLabelString("Ch11 Seg. Length (0 to disable write):");
        deviceField14.setOffsetNid(106);
        jPanel13.add(deviceField14);

        jPanel2.add(jPanel13);

        deviceField15.setIdentifier("");
        deviceField15.setLabelString("Ch12 Seg. Length (0 to disable write):");
        deviceField15.setOffsetNid(114);
        jPanel14.add(deviceField15);

        jPanel2.add(jPanel14);

        deviceField16.setIdentifier("");
        deviceField16.setLabelString("Ch13 Seg. Length (0 to disable write):");
        deviceField16.setOffsetNid(122);
        jPanel15.add(deviceField16);

        jPanel2.add(jPanel15);

        deviceField17.setIdentifier("");
        deviceField17.setLabelString("Ch14 Seg. Length (0 to disable write):");
        deviceField17.setOffsetNid(130);
        jPanel16.add(deviceField17);

        jPanel2.add(jPanel16);

        deviceField18.setIdentifier("");
        deviceField18.setLabelString("Ch15 Seg. Length (0 to disable write):");
        deviceField18.setOffsetNid(138);
        jPanel17.add(deviceField18);

        jPanel2.add(jPanel17);

        deviceField19.setIdentifier("");
        deviceField19.setLabelString("Ch16 Seg. Length (0 to disable write):");
        deviceField19.setOffsetNid(146);
        jPanel18.add(deviceField19);

        jPanel2.add(jPanel18);

        jScrollPane2.setViewportView(jPanel2);

        getContentPane().add(jScrollPane2, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private DeviceButtons deviceButtons1;
    private DeviceField deviceField1;
    private DeviceField deviceField10;
    private DeviceField deviceField11;
    private DeviceField deviceField12;
    private DeviceField deviceField13;
    private DeviceField deviceField14;
    private DeviceField deviceField15;
    private DeviceField deviceField16;
    private DeviceField deviceField17;
    private DeviceField deviceField18;
    private DeviceField deviceField19;
    private DeviceField deviceField2;
    private DeviceField deviceField3;
    private DeviceField deviceField4;
    private DeviceField deviceField5;
    private DeviceField deviceField6;
    private DeviceField deviceField7;
    private DeviceField deviceField8;
    private DeviceField deviceField9;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel10;
    private javax.swing.JPanel jPanel11;
    private javax.swing.JPanel jPanel12;
    private javax.swing.JPanel jPanel13;
    private javax.swing.JPanel jPanel14;
    private javax.swing.JPanel jPanel15;
    private javax.swing.JPanel jPanel16;
    private javax.swing.JPanel jPanel17;
    private javax.swing.JPanel jPanel18;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JPanel jPanel6;
    private javax.swing.JPanel jPanel7;
    private javax.swing.JPanel jPanel8;
    private javax.swing.JPanel jPanel9;
    private javax.swing.JScrollPane jScrollPane2;
    // End of variables declaration//GEN-END:variables
}
