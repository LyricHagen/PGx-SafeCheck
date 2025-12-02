from fastapi import FastAPI, UploadFile, File, Form
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional
import uvicorn

from .pgx_logic import parse_23andme_raw, parse_vcf, analyze

app = FastAPI(
    title="PGx SafeCheck (Demo)",
    description=(
        "Educational PGx demo. Upload a 23andMe-style file or a simple VCF with rsIDs. "
        "Outputs gene–drug association summaries with approximate interpretations.\n\n"
        "⚠️ Not medical advice. Do not change medications based on this demo."
    ),
    version="0.1.0",
)

# Enable CORS for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # or restrict to your frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve the static frontend
app.mount("/", StaticFiles(directory="frontend", html=True), name="static")

@app.post("/api/analyze")
async def api_analyze(
    file: UploadFile = File(...),
    file_format: str = Form(...),  # '23andme' or 'vcf'
):
    text = (await file.read()).decode(errors="ignore")
    try:
        if file_format.lower() == "23andme":
            rs_map = parse_23andme_raw(text)
        elif file_format.lower() == "vcf":
            rs_map = parse_vcf(text)
        else:
            return JSONResponse({"error": "Unsupported file_format; use '23andme' or 'vcf'."}, status_code=400)
        report = analyze(rs_map)
        return JSONResponse(report)
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=500)

if __name__ == "__main__":
    uvicorn.run("backend.app:app", host="0.0.0.0", port=8000, reload=True)
