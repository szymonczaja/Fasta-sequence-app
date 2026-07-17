import streamlit as st
from io import StringIO
from pathlib import Path
from main import FastaParser, SequenceAnalyzer

st.set_page_config(page_title="Sequence Analyzer", layout="wide")
st.markdown("""
<style>
body {background-color: #f5f5f5;}
.block {padding: 20px; background: white; border-radius: 10px; margin-bottom: 20px;}
h2 {color: #4a7a8c;}
</style>
""", unsafe_allow_html=True)

st.title("🧬 Sequence Analyzer")
st.write("Elegancka aplikacja do analizy sekwencji DNA/RNA z Twoimi klasami.")

st.header("📁 Wgraj plik FASTA")
uploaded = st.file_uploader("Wybierz plik FASTA", type=["fasta", "fa", "txt"])

if uploaded:
    temp_path = Path("temp_uploaded.fasta")
    temp_path.write_text(uploaded.getvalue().decode())

    parser = FastaParser(str(temp_path))
    sequences = parser.parse()
    analyzer = SequenceAnalyzer(sequences)

    st.success(f"Wczytano {len(sequences)} sekwencji.")
    st.subheader("🔍 Wczytane sekwencje")
    for seq in sequences:
        with st.expander(seq.header):
            st.code(seq.seq)

    st.header("📊 GC-content")
    gc = analyzer.gc_content()
    for header, pct in gc.items():
        st.write(f"**{header}** → {pct}")

    st.header("🧪 Translacja RNA")
    translations = analyzer.translate_all()
    if translations:
        for header, protein in translations.items():
            st.write(f"**{header}** → {protein}")
    else:
        st.info("Brak sekwencji RNA do translacji.")

    st.header("🔗 Alignment dwóch sekwencji")
    seq_options = [seq.header for seq in sequences]

    col1, col2 = st.columns(2)
    with col1:
        seq1_header = st.selectbox("Sekwencja 1", seq_options)
    with col2:
        seq2_header = st.selectbox("Sekwencja 2", seq_options)

    seq1 = next(s.seq for s in sequences if s.header == seq1_header)
    seq2 = next(s.seq for s in sequences if s.header == seq2_header)

    if st.button("Wykonaj alignment"):
        aligned1, aligned2, score = analyzer.align(seq1, seq2)
        st.subheader("Wynik alignmentu")
        st.code(aligned1)
        st.code(aligned2)
        st.write(f"**Score:** {score}")

    st.header("🧬 Wykrywanie mutacji")
    if st.button("Znajdź mutacje"):
        try:
            muts = analyzer.find_mutations(seq1, seq2)
            if muts:
                for idx, change in muts:
                    st.write(f"Pozycja {idx}: {change}")
            else:
                st.success("Brak mutacji.")
        except Exception as e:
            st.error(str(e))

else:
    st.info("Wgraj plik FASTA, aby rozpocząć analizę.")
